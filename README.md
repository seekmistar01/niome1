# NIOME : Bittensor Subnet(SN55) for Privacy-Safe Genomic Inteligence

Welcome to the NIOME Subnet! This repository contains all the necessary information to get started, understand our subnet architecture, and contribute.

![niome logo image](docs/logo.png)


## Overview

**NIOME** leverages advanced AI models within the Bittensor network to generate large-scale synthetic genomic data that is statistically indistinguishable from real human DNA. Built on a decentralized, incentive-driven framework, NIOME coordinates miners and validators to continuously produce, evaluate, and refine privacy-preserving genomic profiles—enabling reliable, compliant, and scalable genomic research without exposure to real patient data.

## Purpose

Traditionally, genomic research and precision medicine rely on access to real human DNA data. While such data enables highly accurate population and disease modeling, it is heavily constrained by consent requirements, privacy regulations, and the risk of data breaches. As a result, access to large-scale genomic datasets is limited, expensive, and slow—significantly restricting the pace of biomedical research and innovation.

Recent advances in machine learning have enabled the generation of high-fidelity synthetic genomic data that preserves the statistical structure, correlations, and population-level signals of real DNA without containing any identifiable individual information. These data-driven approaches can provide genomic datasets at unlimited scale, with zero privacy risk, and at a fraction of the cost and time required to collect real-world genomic samples.

The NIOME subnet incentivizes the development of novel and robust AI architectures for synthetic genomic data generation and validation. Through the continuous evolution of this subnet, miners are progressively challenged with increasingly complex genomic simulations, enabling higher biological realism, stronger statistical guarantees, and broader applicability across precision medicine and biomedical research.

**System Flow**
1. **Task Generation (Backend → Validators)** The NIOME backend continuously produces genomic simulation tasks. These tasks represent environmental conditions, population parameters, or biological constraints that miners must use to generate synthetic genomic profiles. Validators fetch these tasks directly from the backend. 
2. **Task Distribution (Validators → Miners)** After receiving a task, validators broadcast it to miners on the subnet. Each miner receives the same challenge, ensuring fair and comparable evaluation across participants. 
3. **Synthetic Genome Generation (Miners → Validators)** Miners process the task using their generative models and produce a synthetic genome file. This output must reflect realistic allele frequencies, linkage disequilibrium patterns, and pharmacogenomic variants while preserving privacy. Miners then return their generated genome profile to the validator that issued the task. 
4. **Evaluation and Scoring (Validators → Backend)** Validators evaluate miner submissions using held‑out datasets, statistical fidelity checks, and biological plausibility metrics. Based on performance, validators assign scores that determine miner emissions.

**Features:**

Privacy-Preserving Genomic Generation: NIOME produces high-fidelity synthetic genomic data that preserves real-world statistical structure and biological correlations—without containing any identifiable individual DNA.

Model Evolution: The subnet continuously integrates advances in population genetics, pharmacogenomics, and generative AI to improve realism, coverage, and biological plausibility of synthetic genomes.

Scalable Data Access: Synthetic datasets can be generated at arbitrary scale, enabling large cohort simulations that would be infeasible or unethical using real patient data.

**Core Components:**

- **Miners:** Responsible for running genomic simulation and generative models to produce synthetic DNA profiles, including pharmacogenomically relevant variants such as CYP2D6.

   Incentivized Innovation: Miners are rewarded in $TAO based on the quality, novelty, and statistical fidelity of their generated genomic data. No specialized hardware is required—any participant capable of running a Bittensor miner can contribute.

   Biological Realism: Miner outputs must capture population-level allele frequencies, linkage patterns, and gene–drug response variability observed in real-world datasets.

- **Research Integration:** We systematically update our detection models and methodologies in response to emerging academic research. 

- **Validators:** Responsible for challenging miners with a subsets of environmental data and evaluating miner performance on heldout data.

- **Resource Expansion:** We continuously add new enviromental challenges and data modalities to our subnet in order to evolve our subnet and solve a multitude of distinct problems.

## Guide for Miners and Validators
- [Miner Setup](docs/miner_guide.md) 
- [Validator Setup](docs/validator_guide.md)

## Community
For real-time discussions, community support, and regular updates, <a href="https://discord.com/invite/bittensor">join the bittensor discord</a>. Connect with developers, researchers, and users to get the most out of the NIOME Subnet.

## License
This repository is licensed under the MIT License.
```text
# The MIT License (MIT)
# Copyright © 2024 Opentensor Foundation

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
```