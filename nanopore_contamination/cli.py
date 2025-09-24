"""Command-line interface for the Nanopore contamination pipeline."""
from __future__ import annotations

from typing import Iterable, Optional

from .pipeline import PipelineConfig, build_parser, run_pipeline


def parse_args(argv: Optional[Iterable[str]] = None) -> PipelineConfig:
    """Parse command-line arguments into a :class:`PipelineConfig`."""
    parser = build_parser()
    args = parser.parse_args(argv)
    return PipelineConfig(
        input_fasta=args.input,
        output_dir=args.output,
        search_target=args.search_target,
        force_reannotate=args.force_reannotate,
        email=args.email,
        api_key=args.api_key,
        min_length=args.min_length,
        expect=args.expect,
        hitlist_size=args.hitlist_size,
        delay=args.delay,
        min_identity=args.min_identity,
        min_coverage=args.min_coverage,
        min_matchlen=args.min_matchlen,
        blast_mode=args.blast_mode,
        blastn=args.blastn,
        db=args.db,
        max_target_seqs=args.max_target_seqs,
    )


def main(argv: Optional[Iterable[str]] = None) -> None:
    """Entry point for the console script."""
    config = parse_args(argv)
    run_pipeline(config)


__all__ = ["main", "parse_args"]
