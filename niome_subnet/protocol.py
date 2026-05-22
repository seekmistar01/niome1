"""
Protocol definitions for the Drug Response Prediction Subnet.

This module defines the communication protocols between validators and miners
for drug response prediction tasks using synthetic genomic data.
"""

from typing import Dict, Any, Optional

import bittensor as bt

from niome_subnet.genomics.model import Task

class GenomicsTaskSynapse(bt.Synapse):
    """Protocol for genomics simulation tasks."""

    # Input fields
    task: Optional[Task] = None
    timeout: Optional[float] = None  # Timeout window for submission

    # Output fields
    vcf_content: Optional[str] = None
    cftr_annotations: Optional[Dict[str, Any]] = None
    elapsed_time: Optional[float] = None
    signature: Optional[str] = None  # Cryptographic signature

    def deserialize(self) -> bt.Synapse:
        """Deserialize the GenomicsTaskSynapse Object."""
        return self