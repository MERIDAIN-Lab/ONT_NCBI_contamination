# ONT Adapters BLAST Annotation Pipeline

This repository provides a Python script for identifying contamination in ONT adapter sequences using BLAST searches and automated GenBank metadata retrieval. The pipeline performs the following steps:

1. **Short BLASTn Search** – Runs BLASTn queries for each adapter sequence in `ONT_and_Porechop_adapters_and_barcodes_CLEANED.fasta` and stores the top hits in individual TSV tables.
2. **GenBank Annotation** – Retrieves metadata for each BLAST hit through the NCBI Entrez API and compiles a comprehensive table of annotations.
3. **Summary Reporting** – Generates an Excel workbook combining all annotations and a presence matrix for quick reference. Failed accession queries are listed separately.

The main script is [`blast_annotation_pipeline.py`](blast_annotation_pipeline.py). Results are written to `Combined_BLAST_Report.xlsx` and supporting files in `shortblast_outputs/`.

The workflow was designed for research use and may require adjustment of file paths or Entrez credentials for other environments.
