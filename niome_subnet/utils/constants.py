# ---- General Constants -----
TESTNET_UID = 289
MAINNET_UID = 55

FORWARD_TIMEOUT = 60
MINER_QUERY_K = 5


# ---- Scoring Constants -----
TOP_MINER_COUNT = 10
SCORE_DISTRIBUTION = [0.2, 0.2, 0.2, 0.1, 0.1, 0.05, 0.05, 0.05, 0.025, 0.025]


# ---- Backend Request -----
BASE_URL = "https://niome-api.genomes.io"
MINER_SCORE_URL = f"{BASE_URL}/api/miner_scores"
TASK_URL = f"{BASE_URL}/api/tasks"
GROUND_TRUTH_URL = f"{BASE_URL}/api/tasks/ground_truth"


# ---- Timeout Values -----
TASK_REQUEST_TIMEOUT = 60  # seconds
BASE_DELAY_SECONDS = 2  # seconds
SUBMIT_REQUEST_TIMEOUT = 30  # seconds


# ---- Other Constants -----
MAX_TASK_RETRIES = 3
MAX_SUBMIT_RETRIES = 3

WANDB_MAX_LOGS = 60_000

SCORING_SYSTEM = "top"  # "linear", "top"
BURNING_RATE = 0.02
OWNER_HOTKEY = "5DJ5fT174AY8GzbYHnamYQCJd4cTcj2Zf7ogUvBhry1KfYVd"

BASE_BLOCK_NUMBER = 8034230
INTERVAL_BLOCKS = 2400
FETCHING_BLOCK = 450
VALIDATION_BLOCK = 1050
WEIGHT_SET_BLOCK = 1100
