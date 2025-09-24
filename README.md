# Oxford Nanopore Adapter Contamination Pipeline

This repository accompanies the manuscript describing Oxford Nanopore Technologies (ONT) adapter contamination observed in
GenBank records. It provides a reproducible pipeline for short BLASTn searches, GenBank metadata harvesting, and generation of
summary tables that mirror the analyses in the paper.

All database queries were completed on 23 September 2025.

## Repository contents

- `ONT_and_Porechop_adapters_and_barcodes_CLEANED.fasta` — curated adapter and barcode reference sequences analysed in the study.
- `nanopore_contamination/` — Python package implementing the full annotation workflow.
- `nanopore_adapter_annotation.py` — convenience wrapper for launching the pipeline from the command line.
- `shortblast_outputs/`, `Combined_BLAST_Report.xlsx`, and `partial_annotation.pkl` — published outputs used throughout the
  manuscript.

## Installation

Create a fresh Python environment (Python 3.9 or later) and install the dependencies either with `pip` or `conda`:

```bash
# Using pip
python -m venv .venv
source .venv/bin/activate
pip install biopython pandas python-dateutil openpyxl

# Using conda
conda create -n ont_contamination python=3.11 biopython pandas python-dateutil openpyxl
conda activate ont_contamination
```

If you intend to run BLAST locally, install the NCBI BLAST+ command-line suite and download the relevant databases (e.g. `nt`).

## Usage

1. Provide your NCBI contact details:
   ```bash
   export NCBI_EMAIL="your.email@example.com"
   # Optional: export NCBI_API_KEY="your_ncbi_api_key"
   ```

2. Execute the pipeline:
   ```bash
   python nanopore_adapter_annotation.py --input ONT_and_Porechop_adapters_and_barcodes_CLEANED.fasta \
       --output shortblast_outputs
   ```

   Optional arguments allow you to force re-annotation, direct BLAST to a local database, or run local pairwise searches against
   a custom FASTA. Run `python nanopore_adapter_annotation.py --help` for the full interface.

The command above reproduces the published workflow without modifying the existing outputs. To refresh annotations from GenBank,
pass `--force-reannotate` (note that repeated requests should respect NCBI rate limits).

## Reproducing the manuscript results

1. Clone this repository and install dependencies as described above.
2. Ensure that the BLAST+ executables are discoverable in your `PATH` (or supply `--blastn` and `--db`).
3. Run the pipeline with the provided FASTA file. The script will reuse the packaged `shortblast_outputs/` tables and the
   `partial_annotation.pkl` cache unless `--force-reannotate` is specified.
4. The final Excel workbook (`Combined_BLAST_Report.xlsx`) and intermediate tables match those referenced in the manuscript.

## Notes on reproducibility

- Entrez requests automatically throttle according to the presence of an API key. Adjust `--delay` if you encounter rate-limit
  warnings.
- Local search mode (`--search_target`) requires Biopython's optional `pairwise2` module.
- The pipeline performs no stochastic operations; repeated runs with the same inputs yield identical outputs.
- Sensitive credentials such as API keys should be supplied via environment variables or CLI arguments; no secrets are stored in
  the repository.

For questions related to the manuscript analyses or pipeline usage, please open an issue on the project page.
