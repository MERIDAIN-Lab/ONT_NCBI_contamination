"""Nanopore adapter contamination analysis package."""

from .cli import main
from .pipeline import PipelineConfig, run_pipeline

__all__ = ["main", "PipelineConfig", "run_pipeline"]
