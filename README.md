# Systematic Detection of Oxford Nanopore Adapter and Barcode Contamination in GenBank

## Overview

This repository supports the study **“Systematic Detection of Oxford Nanopore Adapter and Barcode Contamination in GenBank.”**

The objective of this work is to identify residual Oxford Nanopore Technologies (ONT) adapter and barcode sequences embedded in publicly available nucleotide records and to characterize their distribution across taxa and record types.

Incomplete removal of synthetic ONT sequences during library preparation or post-processing can introduce non-biological fragments into genome assemblies. Once deposited, these fragments may be misinterpreted as genuine sequence and can affect downstream analyses, including genome annotation, comparative genomics, primer design, and pathogen surveillance.

This study demonstrates that ONT-derived sequences are recurrent across GenBank records and exhibit structured, non-random patterns consistent with incomplete trimming rather than incidental similarity.

## Key findings

- ONT adapter and barcode fragments were detected across bacterial, viral, and eukaryotic GenBank accessions.
- Rapid-ligation and ligation-based adapters were the most frequently observed sequence classes.
- Multiple ONT-derived sequences frequently co-occurred within the same accession, indicating structured contamination patterns.
- The ONT sequence catalogue produced higher match signal than length- and GC-matched randomized controls, supporting a non-random origin.

## Repository structure

- `scripts/`  
  Core analysis scripts, including the BLAST + Entrez annotation pipeline and control-generation utilities.

- `Sequences_catalogue_and_controls/`  
  Curated ONT adapter and barcode catalogue derived from ONT documentation and Porechop references, along with control FASTA files.

- `ONT_catalogue_shearch_results_March_2026/`  
  Output tables generated from ONT catalogue queries (per-query TSVs and annotated derivatives).

- `Random_gc_matched_controls_search_results_March_2026/`  
  Output tables generated from GC-matched randomized control sequences.

## Methods overview

### Core analysis pipeline

The main workflow is implemented in:

- `scripts/BLAST_pipeline_nanopore_contamination_2.py`

This script performs three main steps:

1. **Sequence similarity search (BLASTn-short)**  
   Each adapter or barcode sequence is queried against a nucleotide database using `blastn-short`, either locally or via remote NCBI BLAST.

2. **Metadata retrieval via Entrez**  
   For each matched accession, metadata are retrieved using the NCBI Entrez API (`Bio.Entrez`).  
   This includes:
   - sequence definition  
   - taxonomy  
   - accession length  
   - molecule type and topology  
   - release year  
   - sequencing technology indicators extracted from GenBank comments  

   The implementation uses XML parsing with fallback to GenBank flatfile parsing to ensure robustness. :contentReference[oaicite:0]{index=0}

3. **Aggregation and reporting**  
   Results are consolidated into a multi-sheet Excel report including:
   - full annotation table  
   - filtered hits  
   - presence matrix across queries  
   - accessions not successfully retrieved  

The pipeline is designed to be resumable through checkpoint files and supports re-annotation when needed.

### Running the analysis

#### Remote BLAST (recommended for reproducibility)

```bash
python scripts/BLAST_pipeline_nanopore_contamination_2.py \
  -i Sequences_catalogue_and_controls/ONT_and_Porechop_adapters_and_barcodes_CLEANED.fasta \
  -o shortblast_outputs \
  --blast-mode remote \
  --db nt \
  --email your.email@institution.edu
