from pydantic import BaseModel
from typing import Any, Dict, Generic, Optional, TypeVar


class TaskInput(BaseModel):
    read1_fastq: str
    read2_fastq: str


class TaskOutputSpec(BaseModel):
    format: str
    required_fields: list[str]


class TaskGenomeContext(BaseModel):
    chromosome: str
    region: str
    gene: str


class Task(BaseModel):
    task_id: str
    version: str
    type: str

    input: TaskInput
    output_spec: TaskOutputSpec
    genome_context: TaskGenomeContext
    expected_variant_count: int


class GroundTruth(BaseModel):
    truth_vcf: str
    ref: str
    cftr2_annotations: str


class MinerSubmission(BaseModel):
    uid: int
    vcf_content: str
    response_time: float
    cftr_annotations: Optional[Dict[str, Any]] = None


class ValidationContext:
    """Container for all validation metadata."""
    miner_uid: int
    miner_hotkey: str
    validator_uid: int
    validator_hotkey: str

    def __init__(self, miner_uid: int = 0, miner_hotkey: str = "", validator_uid: int = 0, validator_hotkey: str = ""):
        self.miner_uid = miner_uid
        self.miner_hotkey = miner_hotkey
        self.validator_uid = validator_uid
        self.validator_hotkey = validator_hotkey


class TaskPayload(BaseModel):
    """Payload structure for task generation requests."""
    timestamp: float
    hotkey: str
    uuid: str
    netuid: str


PayloadType = TypeVar('PayloadType', bound=BaseModel)

class SignedRequest(BaseModel, Generic[PayloadType]):
    """Generic signed request structure."""
    payload: PayloadType
    payload_raw: str
    signature: str


class MinerScore(BaseModel):
    uid: int
    precision: float
    recall: float
    f1_score: float
    response_time: float
    vcf_score: float
    annotation_score: float
    final_score: float
    log: str


class MinerScoreDto(BaseModel):
    """Data transfer object for miner score submission."""
    task_id: str
    uid: int
    hotkey: str
    precision: float
    recall: float
    f1_score: float
    response_time: float
    vcf_score: float
    annotation_score: float
    final_score: float
    log: str
    weight: float
