# The MIT License (MIT)
# Copyright © 2023 Yuma Rao
# Copyright © 2025 Genomes.io

# Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated
# documentation files (the “Software”), to deal in the Software without restriction, including without limitation
# the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software,
# and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

# The above copyright notice and this permission notice shall be included in all copies or substantial portions of
# the Software.

# THE SOFTWARE IS PROVIDED “AS IS”, WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO
# THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
# THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
# OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
# DEALINGS IN THE SOFTWARE.

import aiohttp
import asyncio
import bittensor as bt
import copy
import json
import numpy as np
import niome_subnet.utils.constants as config
import os
import time
import urllib.request

from typing import Optional
from niome_subnet.genomics.model import GroundTruth, Task, MinerSubmission
from niome_subnet.genomics.scoring import create_mapping_file, score
from niome_subnet.protocol import GenomicsTaskSynapse
from niome_subnet.utils import get_miner_uids

from niome_subnet.utils.constants import BASE_BLOCK_NUMBER, BURNING_RATE, FETCHING_BLOCK, INTERVAL_BLOCKS, VALIDATION_BLOCK

sem = asyncio.Semaphore(config.MINER_QUERY_K)


async def fetch_task(self) -> Task:
    """Generate a synthetic genomic simulation task with retry logic and fallback."""
    payload = {}
    timestamp = str(time.time())
    canonical = json.dumps({
        'payload': '{}',
        'hotkey': self.wallet.hotkey.ss58_address,
        'netuid': str(self.netuid),
        'timestamp': timestamp,
    }, separators=(',', ':'), sort_keys=True)

    signature = self.wallet.hotkey.sign(canonical).hex()

    for attempt in range(1, config.MAX_TASK_RETRIES + 1):
        try:
            async with aiohttp.ClientSession() as client:
                async with client.post(
                    config.TASK_URL,
                    headers=self.build_signature_headers(
                        signature=signature,
                        hotkey=self.wallet.hotkey.ss58_address,
                        timestamp=timestamp,
                        netuid=str(self.netuid),
                    ),
                    json=payload,
                    timeout=aiohttp.ClientTimeout(total=config.TASK_REQUEST_TIMEOUT),
                ) as response:
                    if response.status != 201:
                        raise RuntimeError(
                            f"Backend returned status {response.status}"
                        )

                    data = await response.json()
                    
                    task_url = data.get("task_url", "")

                    if not task_url:
                        raise RuntimeError("Invalid response from backend")

                    task = await fetch_task_by_url(task_url)

                    return task
        except Exception as e:
            bt.logging.error(f"Error on generating task (attempt {attempt}): {e}")
            if attempt < config.MAX_TASK_RETRIES:
                delay = config.BASE_DELAY_SECONDS * (2 ** (attempt - 1))
                bt.logging.info(f"Retrying in {delay} seconds...")
                await asyncio.sleep(delay)
            else:
                bt.logging.error("All retries failed, returning fallback sample data")
                raise e

async def fetch_ground_truth(self) -> GroundTruth:
    """Generate a synthetic genomic simulation task with retry logic and fallback."""
    payload = {}
    timestamp = str(time.time())
    canonical = json.dumps({
        'payload': '{}',
        'hotkey': self.wallet.hotkey.ss58_address,
        'netuid': str(self.netuid),
        'timestamp': timestamp,
    }, separators=(',', ':'), sort_keys=True)

    signature = self.wallet.hotkey.sign(canonical).hex()

    for attempt in range(1, config.MAX_TASK_RETRIES + 1):
        try:
            async with aiohttp.ClientSession() as client:
                async with client.post(
                    config.GROUND_TRUTH_URL,
                    headers=self.build_signature_headers(
                        signature=signature,
                        hotkey=self.wallet.hotkey.ss58_address,
                        timestamp=timestamp,
                        netuid=str(self.netuid),
                    ),
                    json=payload,
                    timeout=aiohttp.ClientTimeout(total=config.TASK_REQUEST_TIMEOUT),
                ) as response:
                    if response.status != 201:
                        raise RuntimeError(
                            f"Backend returned status {response.status}"
                        )

                    data = await response.json()
                    
                    ground_truth_url = data.get("ground_truth_url", "")

                    if not ground_truth_url:
                        raise RuntimeError("Invalid response from backend")

                    ground_truth = await fetch_ground_truth_by_url(ground_truth_url)

                    return ground_truth
        except Exception as e:
            bt.logging.error(f"Error on fetching ground truth (attempt {attempt}): {e}")
            if attempt < config.MAX_TASK_RETRIES:
                delay = config.BASE_DELAY_SECONDS * (2 ** (attempt - 1))
                bt.logging.info(f"Retrying in {delay} seconds...")
                await asyncio.sleep(delay)
            else:
                bt.logging.error("All retries failed, returning fallback sample data")
                raise e

async def fetch_task_by_url(task_url: str) -> Task:
    """Fetch task details from the given URL."""
    try:
        def _fetch():
            with urllib.request.urlopen(task_url, timeout=config.TASK_REQUEST_TIMEOUT) as resp:
                return json.loads(resp.read())
        data = await asyncio.to_thread(_fetch)
        return Task(**data)
    except Exception as e:
        bt.logging.error(f"Error fetching task details from {task_url}: {type(e).__name__}: {e}")
        raise

async def fetch_ground_truth_by_url(ground_truth_url: str) -> GroundTruth:
    """Fetch ground truth data from the given URL."""
    try:
        def _fetch():
            with urllib.request.urlopen(ground_truth_url, timeout=config.TASK_REQUEST_TIMEOUT) as resp:
                return json.loads(resp.read())
        data = await asyncio.to_thread(_fetch)
        return GroundTruth(**data)
    except Exception as e:
        bt.logging.error(f"Error fetching ground truth from {ground_truth_url}: {type(e).__name__}: {e}")
        raise

async def query_axon(self, axon, synapse) -> Optional[GenomicsTaskSynapse]:
    """Query a single axon and return the response."""
    try:
        start_time = time.perf_counter()
        response: GenomicsTaskSynapse = await self.dendrite.forward(
            axons=axon, synapse=synapse, deserialize=False, timeout=config.FORWARD_TIMEOUT
        )
        if response is not None:
            response.elapsed_time = time.perf_counter() - start_time
            return response
    except Exception as e:
        bt.logging.error(f"Error querying axon {axon}: {e}")
        return None

async def collect_miners_responses(self):
    bt.logging.info("Collecting miners' responses...")
    try:
        os.makedirs("data", exist_ok=True)
        os.makedirs("vcfs", exist_ok=True)
        miner_uids = get_miner_uids(self)
        np.random.shuffle(miner_uids)

        miner_task = await fetch_task(self)
        bt.logging.info(f"Fetched task: {miner_task.model_dump()}")
        task = copy.deepcopy(miner_task)
        self.task_id = task.task_id

        # Download task reads
        urllib.request.urlretrieve(task.input.read1_fastq, "data/read_1.fq")
        task.input.read1_fastq = "data/read_1.fq"
        urllib.request.urlretrieve(task.input.read2_fastq, "data/read_2.fq")
        task.input.read2_fastq = "data/read_2.fq"

        for uid in miner_uids:
            synapse = GenomicsTaskSynapse(task=miner_task, timeout=config.FORWARD_TIMEOUT)

            axon = self.metagraph.axons[uid]
            if axon.ip == '0.0.0.0':
                continue

            response = await query_axon(self, axon, synapse)
            if response is None or response.vcf_content is None:
                continue

            lines = response.vcf_content.splitlines()
            variant_count = 0
            for line in lines:
                if not line.startswith("#"):
                    variant_count += 1

            with open(f"vcfs/{uid}.vcf", "w") as f:
                vcf_content = f"##response_time={response.elapsed_time}\n" + response.vcf_content
                f.write(vcf_content)

            if response.cftr_annotations is not None:
                with open(f"vcfs/{uid}.annotations.json", "w") as f:
                    json.dump(response.cftr_annotations, f)

        self.is_validating = False
        bt.logging.info("Finished collecting responses.")
    except Exception as e:
        bt.logging.error(f"Error during fetching process: {e}")
    finally:
        self.is_fetching = False

async def run_validation(self):
    try:
        bt.logging.info("Validating miners' responses...")

        ground_truth = await fetch_ground_truth(self)
        bt.logging.info(f"Fetched ground truth")

        # Download ground truth data first (ref needed by create_mapping_file)
        urllib.request.urlretrieve(ground_truth.truth_vcf, "data/truth.vcf")
        ground_truth.truth_vcf = "data/truth.vcf"
        urllib.request.urlretrieve(ground_truth.ref, "data/ref.fa")
        ground_truth.ref = "data/ref.fa"
        urllib.request.urlretrieve(ground_truth.cftr2_annotations, "data/cftr2_annotations.json")
        ground_truth.cftr2_annotations = "data/cftr2_annotations.json"

        bam = create_mapping_file(ground_truth.ref, "data/read_1.fq", "data/read_2.fq")

        final_scores: list[MinerSubmission] = []

        for vcf_file in os.listdir("vcfs"):
            if not vcf_file.endswith(".vcf"):
                continue
            uid = int(os.path.splitext(vcf_file)[0])
            with open(f"vcfs/{vcf_file}") as f:
                lines = f.readlines()
            response_time = None
            vcf_lines = []
            for line in lines:
                if line.startswith("##response_time="):
                    response_time = float(line.strip().split("=", 1)[1])
                else:
                    vcf_lines.append(line)
            vcf_content = "".join(vcf_lines)

            cftr_annotations = None
            annotations_path = f"vcfs/{os.path.splitext(vcf_file)[0]}.annotations.json"
            if os.path.exists(annotations_path):
                with open(annotations_path) as f:
                    cftr_annotations = json.load(f)

            miner_score = score(
                MinerSubmission(
                    uid=uid,
                    vcf_content=vcf_content,
                    response_time=response_time,
                    cftr_annotations=cftr_annotations,
                ), ground_truth, bam)

            final_scores.append(miner_score)

        bt.logging.info(f"Scores: {[(score.uid, score.vcf_score, score.annotation_score, score.final_score) for score in final_scores]}")
        self.set_weights(final_scores, self.task_id)
    except Exception as e:
        bt.logging.error(f"Error validating miners' vcf: {e}")

async def forward(self):
    """
    The forward function is called by the validator every time step.

    It is responsible for querying the network and scoring the responses.

    Args:
        self (:obj:`bittensor.neuron.Neuron`): The neuron object which contains all the necessary state for the validator.

    """
    try:
        if BURNING_RATE == 1.0 and not self.is_validating:
            self.is_validating = True
            self.set_weights([], "")
            self.subtensor.set_weights(
                wallet=self.wallet,
                netuid=self.config.netuid,
                uids=self.uids,
                weights=self.weights,
                wait_for_finalization=False,
                wait_for_inclusion=False,
            )
        else:
            if (self.block - BASE_BLOCK_NUMBER) % INTERVAL_BLOCKS == FETCHING_BLOCK and not self.is_fetching:
                self.is_fetching = True
                self.are_weights_committed = False
                asyncio.create_task(collect_miners_responses(self))
            elif (self.block - BASE_BLOCK_NUMBER) % INTERVAL_BLOCKS == VALIDATION_BLOCK and not self.is_validating:
                self.is_validating = True
                asyncio.create_task(run_validation(self))
    except Exception as e:
        bt.logging.error(f"Error during forward step: {e}")

    time.sleep(5)
