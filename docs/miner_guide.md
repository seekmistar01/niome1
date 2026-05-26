##  Miner Setup Guide

To contribute as a Miner in the NIOME(SN55) and earn **TAO (**$\tau$**)** emissions, you should first set up your machine, install the necessary packages, and register your identity on the Bittensor.

### 1. Prerequisites

Before starting, ensure your system meets the minimum requirements and has the core dependencies installed.

* **Operating System:** Ubuntu 22.04 or similar Linux distribution is generally recommended for optimal compatibility. Mining is not supported on Windows.
* **Python:** Python 3.12
* **Git:** cloning the repository
* **Hardware:** 
   - vCPU : 16
   - GPU : necessary, up to your generation
   - Memory : 16GB minimum
   - 3rd Party API : unnecessary
   - Port Forwarding : standard

### 2. Environment Setup

This section walks you through cloning the subnet-niome repository and installing the required packages.

1. **Clone the Repository:**
   **Bash**

   ```
   git clone https://github.com/genomesio/subnet-niome.git
   cd subnet-niome
   ```
2. **Create a Virtual Environment (Recommended):**
   **Bash**

   ```
   python3 -m venv venv
   source venv/bin/activate
   ```
3. **Install Dependencies:** Install the required Python packages and register the local package for execution.
   **Bash**

   ```
   python3 -m pip install -r requirements.txt
   ```

### 3. Wallet Creation and Registration

You must create a Bittensor wallet to hold your TAO and Alpha tokens, and to register your hotkey with the subnet.

1. **Install Bittensor CLI:**
   **Bash**

   ```
   python3 -m pip install bittensor-cli
   ```
2. **Create a Coldkey (Primary Wallet):** The coldkey is your secure, offline store of funds. Choose a secure name
   **Bash**

   ```
   btcli wallet new_coldkey --wallet.name your_coldkey
   ```
3. **Create a Hotkey (Miner Identity):** The hotkey is used to sign transactions, run the miner, and receive emissions. It is connected to your coldkey.
   **Bash**

   ```
   btcli wallet new_hotkey --wallet.name your_coldkey --wallet.hotkey your_hotkey
   ```
4. **Fund Your Coldkey:** Transfer a small amount of TAO to your coldkey to cover registration fees, which fluctuate based on subnet competition.
5. **Register Your Hotkey to Subnet 55:** Register your hotkey to secure a UID (Unique Identifier) on the NIOME subnet. The Network ID for NIOME is 55.
   **Bash**

   ```
   btcli subnet register --netuid 55 --wallet.name your_coldkey --wallet.hotkey your_hotkey
   ```

### 4. Install miner tools and reference data

From the project root:

```bash
bash scripts/setup_miner_env.sh
bash scripts/verify_miner_env.sh
```

This installs `bwa`, `samtools`, `bcftools`, `tabix`, `screen`, Java, GATK 4.6.2.0, a Python venv with `requirements.txt`, and downloads/indexes `data/chr7.fa` (GRCh38).

### 5. Sensitive multi-caller strategy (recommended)

The miner uses `niome_subnet/strategy_auto.py` by default for NIOME’s low-coverage,
allele-imbalance simulation (15–85% haplotype fractions). It adds bcftools, FreeBayes
(Docker), and a MAPQ-relaxed weak scanner on top of GATK + ClinVar panel probes.

Environment variables (optional):

| Variable | Default | Meaning |
|----------|---------|---------|
| `NIOME_USE_STRATEGY_AUTO` | `1` | Set `0` to disable strategy merge |

Live miners load these from `scripts/miner_env.sh` (sourced by `scripts/start_miners_screen.sh`). Restart after code changes:

```bash
bash scripts/restart_miners.sh
```
| `NIOME_STRATEGY_SENSITIVITY` | `sensitive` | `balanced`, `sensitive`, or `aggressive` |
| `NIOME_STRATEGY_SUBMIT_TOP_N` | `30` | Cap strategy-aligned submission size |
| `NIOME_STRATEGY_SUBMIT_ALLOWLIST` | `1` | Submit from strategy rank order (selected set) |
| `NIOME_STRATEGY_USE_RANK_GT` | `1` | Use GT from `ranked_candidates.tsv` in submission |
| `NIOME_STRATEGY_HYBRID_HOM_GT` | `1` | Upgrade to `1/1` when BAM AF supports hom (validator GT score) |
| `NIOME_STRATEGY_DEEPVARIANT` | `0` | Set `1` to also run DeepVariant in Docker |
| `NIOME_EVIDENCE_MIN_ALT_DEPTH` | `1` | Min alt reads when strategy is on |
| `NIOME_EVIDENCE_MIN_AF` | `0.08` | Min allele fraction for evidence |
| `NIOME_EVIDENCE_SNV_MIN_AF` | `0.10` | SNV AF floor with strategy on |

Strategy artifacts are written under `work/<task>/<instance>/strategy_auto/`.

### 6. Running the Miner

Once your hotkey is registered, you can start your Miner.

**Single miner:**

```bash
export PYTHONPATH="$(pwd)"
export NIOME_GATK="$(pwd)/tools/gatk/gatk"
./venv/bin/python neurons/miner.py \
  --netuid 55 \
  --subtensor.network finney \
  --wallet.name your_coldkey \
  --wallet.hotkey your_hotkey \
  --axon.port 50007 \
  --logging.info
```

**Four miners in screen** (set wallet coldkeys before starting):

```bash
export NIOME_WALLET_MAIN=your_coldkey
export NIOME_WALLET_ALT=your_other_coldkey   # optional, default seekmistar3
bash scripts/start_miners_screen.sh
```

Use **`screen -r niome_m1`** (etc.) to attach to a session.

### 7. Keep it Running

Validators reward only active, responsive miners. Use `screen`, `tmux`, or `pm2` so processes survive disconnects.