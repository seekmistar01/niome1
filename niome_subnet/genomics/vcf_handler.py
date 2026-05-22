import json
import time
import requests
import bittensor as bt
from niome_subnet.genomics.model import MinerScoreDto
from niome_subnet.utils.constants import (
    MINER_SCORE_URL,
    BASE_DELAY_SECONDS,
    MAX_SUBMIT_RETRIES,
)


def submit_validation_result(
    self,
    miner_scores: list[MinerScoreDto],
):
    try:
        for attempt in range(1, MAX_SUBMIT_RETRIES + 1):
            try:
                payload = {
                    "miner_scores": json.dumps([score.model_dump() for score in miner_scores], indent=2)
                }
                timestamp = str(time.time())
                canonical = json.dumps({
                    'payload': json.dumps(payload, separators=(',', ':'), sort_keys=True),
                    'hotkey': self.wallet.hotkey.ss58_address,
                    'netuid': str(self.netuid),
                    'timestamp': timestamp,
                }, separators=(',', ':'), sort_keys=True)

                signature = self.wallet.hotkey.sign(canonical).hex()

                response = requests.post(
                    MINER_SCORE_URL,
                    headers=self.build_signature_headers(
                        signature=signature,
                        hotkey=self.wallet.hotkey.ss58_address,
                        timestamp=timestamp,
                        netuid=str(self.netuid),
                    ),
                    json=payload
                )
                if response.status_code == 200:
                    bt.logging.info("Validation result submitted successfully")
                    break
                else:
                    bt.logging.error(
                        f"Backend submission failed: {response.status_code} | {response.text}"
                    )
                    if attempt == MAX_SUBMIT_RETRIES:
                        bt.logging.error("All retries failed, giving up.")
            except Exception as e:
                bt.logging.error(f"Error submitting validation result (attempt {attempt}): {e}")
                if attempt == MAX_SUBMIT_RETRIES:
                    bt.logging.error("All retries failed, giving up.")
            if attempt < MAX_SUBMIT_RETRIES:
                delay = BASE_DELAY_SECONDS * (2 ** (attempt - 1))
                bt.logging.info(f"Retrying in {delay} seconds...")
                time.sleep(delay)
    except Exception as e:
        bt.logging.error(f"Unexpected error in submit_validation_result: {e}")
        
