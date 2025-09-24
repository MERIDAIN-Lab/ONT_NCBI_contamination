"""Core pipeline for annotating Oxford Nanopore adapter contamination."""
from __future__ import annotations

import argparse
import logging
import os
import re
import socket
import subprocess
import tempfile
import time
from contextlib import closing
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple, Set

import pandas as pd
from dateutil import parser as dateparser
from Bio import Entrez, SeqIO
import xml.etree.ElementTree as ET
import xml.parsers.expat as expat

try:  # Optional import for local alignment mode
    from Bio import pairwise2
    HAVE_PAIRWISE2 = True
except Exception:  # pragma: no cover - optional dependency
    HAVE_PAIRWISE2 = False

TSV_HEADER: List[str] = [
    "Accession", "Length", "E-value", "Definition",
    "Organism", "Taxonomy", "Release Year",
    "query_start", "query_end", "subject_start", "subject_end",
    "Align_Length", "Query_Length", "Identity(%)", "Coverage(%)", "Pass_Filter",
]


@dataclass
class PipelineConfig:
    """Configuration for the adapter contamination pipeline."""

    input_fasta: str = "ONT_and_Porechop_adapters_and_barcodes_CLEANED.fasta"
    output_dir: str = "shortblast_outputs"
    search_target: Optional[str] = None
    force_reannotate: bool = False
    email: str = os.environ.get("NCBI_EMAIL", "")
    api_key: str = os.environ.get("NCBI_API_KEY", "")
    min_length: int = 8
    expect: float = 10.0
    hitlist_size: int = 25
    delay: Optional[float] = None
    min_identity: float = 95.0
    min_coverage: float = 80.0
    min_matchlen: int = 16
    blast_mode: str = "local"
    blastn: str = "blastn"
    db: str = "nt"
    max_target_seqs: int = 100_000_000
    final_excel: str = "Combined_BLAST_Report.xlsx"
    not_found_file: str = "accessions_not_retrieved.txt"
    partial_pickle: str = "partial_annotation.pkl"

    @property
    def effective_delay(self) -> float:
        """Return the throttling delay for Entrez requests."""
        if self.delay is not None:
            return float(self.delay)
        return 0.34 if self.api_key else 0.5


def configure_logging() -> None:
    """Configure module-level logging."""
    logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")


def configure_entrez(config: PipelineConfig) -> None:
    """Populate Entrez configuration from the pipeline settings."""
    Entrez.email = config.email or "your.email@example.com"
    if config.api_key:
        Entrez.api_key = config.api_key


def sanitize_id(identifier: str, maxlen: int = 64) -> str:
    """Return a filesystem-safe identifier."""
    clean = re.sub(r"[^A-Za-z0-9]+", "_", identifier).strip("_")
    if len(clean) > maxlen:
        return clean[:maxlen]
    return clean


def tsv_path_for(seq_id: str, seq: str, output_dir: str) -> str:
    """Generate the per-query TSV path for a sequence."""
    clean_id = sanitize_id(seq_id)
    suffix_core = seq if len(seq) <= 30 else f"{seq[:15]}...{seq[-10:]}"
    suffix = sanitize_id(suffix_core, maxlen=40)
    filename = f"{clean_id}__{suffix}_description_table.tsv"
    return os.path.join(output_dir, filename)


def write_tsv_header(path: str) -> None:
    """Write the TSV header for BLAST hits."""
    with open(path, "w", encoding="utf-8") as handle:
        handle.write("\t".join(TSV_HEADER) + "\n")


def append_tsv_row(path: str, row: Sequence[str]) -> None:
    """Append a row to the TSV file."""
    with open(path, "a", encoding="utf-8") as handle:
        handle.write("\t".join(row) + "\n")


def parse_year_from_xml(root_el: ET.Element) -> Optional[int]:
    """Extract the release year from an Entrez XML record."""
    dt = root_el.findtext(".//GBSeq_update-date") or root_el.findtext(".//GBSeq_create-date")
    if not dt:
        return None
    try:
        return dateparser.parse(dt).year
    except Exception:
        match = re.search(r"\d{2}-[A-Za-z]{3}-(\d{4})", dt)
        if match:
            try:
                return int(match.group(1))
            except ValueError:
                return None
    return None


def parse_taxonomy_from_flatfile(gb_text: str) -> str:
    """Parse the taxonomy string from a GenBank flatfile."""
    if not gb_text:
        return ""
    match = re.search(r"^ORGANISM[ \t]+[^\n]*\n((?:[ \t]{2,}.+\n)+)", gb_text, flags=re.MULTILINE)
    if not match:
        fallback = re.search(r"^TAXONOMY[ \t]+(.+)$", gb_text, flags=re.MULTILINE)
        if fallback:
            return fallback.group(1).strip().rstrip(".")
        return ""
    block = match.group(1)
    lineage = " ".join(line.strip() for line in block.splitlines()).rstrip(".")
    return lineage


def parse_organism_from_flatfile(gb_text: str) -> str:
    """Parse the organism line from a GenBank flatfile."""
    if not gb_text:
        return "NA"
    match = re.search(r"^ORGANISM[ \t]+(.+)$", gb_text, flags=re.MULTILINE)
    return match.group(1).strip() if match else "NA"


def parse_assembly_lines_from_comment(comment_text: str) -> str:
    """Extract assembly-related lines from a GenBank comment block."""
    if not comment_text:
        return ""
    lines = [line.strip() for line in comment_text.splitlines() if "assembly" in line.lower()]
    return "; ".join(lines)


def detect_sequencing_technology(comment_text: str) -> str:
    """Identify Oxford Nanopore sequencing mentions in the comment."""
    if not comment_text:
        return ""
    lowered = comment_text.lower()
    nanopore_keywords = [
        "nanopore", "ont", "minion", "gridion", "promethion", "flongle", "minknow", "oxford nanopore",
    ]
    if any(keyword in lowered for keyword in nanopore_keywords):
        return "Oxford Nanopore"
    return ""


def entrez_fetch_xml(acc: str, retries: int = 3) -> Optional[str]:
    """Fetch an accession record in XML format from Entrez."""
    for attempt in range(1, retries + 1):
        try:
            with closing(Entrez.efetch(db="nuccore", id=acc, retmode="xml")) as handle:
                return handle.read()
        except Exception as error:  # pragma: no cover - network failure handling
            logging.warning("efetch XML failed for %s (attempt %s): %s", acc, attempt, error)
            time.sleep(1.5 * attempt)
    return None


def entrez_fetch_flat(acc: str, retries: int = 2) -> Optional[str]:
    """Fetch an accession record in GenBank flatfile format from Entrez."""
    for attempt in range(1, retries + 1):
        try:
            with closing(Entrez.efetch(db="nuccore", id=acc, rettype="gb", retmode="text")) as handle:
                return handle.read()
        except Exception as error:  # pragma: no cover - network failure handling
            logging.warning("efetch flatfile failed for %s (attempt %s): %s", acc, attempt, error)
            time.sleep(1.5 * attempt)
    return None


def get_minimal_meta(
    acc: str,
    config: PipelineConfig,
    cache: Dict[str, Dict[str, str]],
) -> Dict[str, str]:
    """Retrieve minimal metadata for an accession using XML with flatfile fallbacks."""
    if acc in cache:
        return cache[acc]

    meta: Dict[str, str] = {"Organism": "NA", "Taxonomy": "", "Release Year": ""}
    xml_data = entrez_fetch_xml(acc, retries=3)
    if xml_data:
        try:
            root = ET.fromstring(xml_data)
            gbseq = root.find(".//GBSeq")
            if gbseq is not None:
                organism = gbseq.findtext("GBSeq_organism") or "NA"
                taxonomy = gbseq.findtext("GBSeq_taxonomy") or ""
                year = parse_year_from_xml(gbseq)
                meta["Organism"] = organism
                meta["Taxonomy"] = taxonomy
                meta["Release Year"] = str(year) if year is not None else ""
        except (ET.ParseError, expat.ExpatError) as error:  # pragma: no cover - XML edge cases
            logging.warning("XML parse error for %s: %s", acc, error)

    if not meta["Taxonomy"] or meta["Organism"] == "NA" or not meta["Release Year"]:
        flat = entrez_fetch_flat(acc, retries=2)
        if flat:
            if not meta["Organism"] or meta["Organism"] == "NA":
                meta["Organism"] = parse_organism_from_flatfile(flat) or "NA"
            if not meta["Taxonomy"]:
                meta["Taxonomy"] = parse_taxonomy_from_flatfile(flat) or ""
            if not meta["Release Year"]:
                match = re.search(r"Date:\s*(\d{2}-[A-Za-z]{3}-\d{4})", flat)
                if match:
                    try:
                        meta["Release Year"] = str(dateparser.parse(match.group(1)).year)
                    except Exception:  # pragma: no cover - unexpected date format
                        pass
    cache[acc] = meta
    time.sleep(config.effective_delay)
    return meta


def run_blast_cli_return_rows(
    seq: str,
    config: PipelineConfig,
    cache: Dict[str, Dict[str, str]],
) -> List[List[str]]:
    """Run the BLAST CLI for a short query sequence and return parsed rows."""
    with tempfile.NamedTemporaryFile(mode="w", suffix=".fa", delete=False) as temp_fasta:
        temp_fasta.write(">q\n")
        temp_fasta.write(seq + "\n")
        query_path = temp_fasta.name

    cmd: List[str] = [
        config.blastn,
        "-task", "blastn-short",
        "-query", query_path,
        "-word_size", "7",
        "-evalue", str(config.expect),
        "-dust", "no",
        "-soft_masking", "false",
        "-db", config.db,
        "-outfmt", "6 sacc evalue pident length qlen qstart qend sstart send qcovs slen",
        "-max_target_seqs", str(config.max_target_seqs),
    ]
    if config.blast_mode == "remote":
        cmd.append("-remote")

    try:
        result = subprocess.run(cmd, check=True, capture_output=True, text=True)
        lines = [line for line in result.stdout.splitlines() if line.strip()]
    except subprocess.CalledProcessError as error:
        logging.error("blastn failed: %s\nSTDERR: %s", error, error.stderr)
        lines = []
    finally:
        try:
            os.unlink(query_path)
        except OSError:
            pass

    rows: List[List[str]] = []
    for line in lines:
        parts = line.split("\t")
        if len(parts) < 11:
            continue
        sacc, evalue, pident, alen, qlen, qstart, qend, sstart, send, qcovs, slen = parts[:11]
        meta = get_minimal_meta(sacc, config, cache)
        organism = meta.get("Organism") or "NA"
        taxonomy = meta.get("Taxonomy", "")
        year = meta.get("Release Year", "")
        try:
            identity_pct = float(pident)
            coverage_pct = float(qcovs)
            align_length = int(float(alen))
        except ValueError:
            identity_pct, coverage_pct, align_length = 0.0, 0.0, 0
        pass_filter = (
            align_length >= int(config.min_matchlen)
            and identity_pct >= float(config.min_identity)
            and coverage_pct >= float(config.min_coverage)
        )
        row = [
            sacc,
            slen,
            evalue,
            "",
            organism,
            taxonomy,
            year,
            qstart,
            qend,
            sstart,
            send,
            alen,
            qlen,
            f"{identity_pct:.2f}",
            f"{coverage_pct:.2f}",
            "YES" if pass_filter else "NO",
        ]
        rows.append(row)
    return rows


def load_target_sequences(path: str) -> Dict[str, str]:
    """Load sequences from a FASTA file for local search mode."""
    return {record.id: str(record.seq) for record in SeqIO.parse(path, "fasta")}


def run_step1(config: PipelineConfig, cache: Dict[str, Dict[str, str]]) -> None:
    """Execute BLAST or local searches and cache minimal metadata per query."""
    print("[STEP 1] Generating TSVs per query (skip existing)")

    target_seqs: Dict[str, str] = {}
    if config.search_target:
        if not HAVE_PAIRWISE2:
            raise RuntimeError("pairwise2 not available. Install Biopython to use --search_target.")
        target_seqs = load_target_sequences(config.search_target)

    for record in SeqIO.parse(config.input_fasta, "fasta"):
        seq_id = record.id
        sequence = str(record.seq)
        if len(sequence) < int(config.min_length):
            print(f"[SKIP] {seq_id}: sequence too short ({len(sequence)} bp)")
            continue

        output_path = tsv_path_for(seq_id, sequence, config.output_dir)
        if os.path.exists(output_path):
            print(f"[SKIP] {seq_id}: already processed -> {os.path.basename(output_path)}")
            continue

        print(f"[INFO] Processing {seq_id} ({len(sequence)} bp) -> {os.path.basename(output_path)}")
        write_tsv_header(output_path)

        if config.search_target:
            best: Optional[Tuple[str, float]] = None
            for target_id, target_seq in target_seqs.items():
                try:
                    alignments = pairwise2.align.localms(sequence, target_seq, 2, -1, -0.5, -0.1)
                except Exception as error:  # pragma: no cover - alignment edge cases
                    logging.warning("Local align failed for %s vs %s: %s", seq_id, target_id, error)
                    alignments = []
                if alignments:
                    score = alignments[0].score
                    if best is None or score > best[1]:
                        best = (target_id, score)
            if best:
                row = [
                    best[0],
                    str(len(target_seqs[best[0]])),
                    "NA",
                    "Local_match",
                    "NA",
                    "",
                    "",
                    "NA",
                    "NA",
                    "NA",
                    "NA",
                    "NA",
                    str(len(sequence)),
                    "NA",
                    "NA",
                    "NA",
                ]
                append_tsv_row(output_path, row)
            else:
                append_tsv_row(
                    output_path,
                    [
                        "NA",
                        "NA",
                        "NA",
                        "No_match_found",
                        "NA",
                        "",
                        "",
                        "NA",
                        "NA",
                        "NA",
                        "NA",
                        "NA",
                        str(len(sequence)),
                        "NA",
                        "NA",
                        "NA",
                    ],
                )
            print(f"[OK] Local search saved: {os.path.basename(output_path)}")
            continue

        rows = run_blast_cli_return_rows(sequence, config, cache)
        if not rows:
            append_tsv_row(
                output_path,
                [
                    "NA",
                    "NA",
                    "NA",
                    "No_hits",
                    "NA",
                    "",
                    "",
                    "NA",
                    "NA",
                    "NA",
                    "NA",
                    "NA",
                    str(len(sequence)),
                    "NA",
                    "NA",
                    "NA",
                ],
            )
        else:
            for row in rows:
                append_tsv_row(output_path, row)
        print(f"[OK] BLAST results saved: {os.path.basename(output_path)}")


def collect_hit_tables(output_dir: str) -> List[pd.DataFrame]:
    """Load all per-query TSV tables from the output directory."""
    dataframes: List[pd.DataFrame] = []
    for filename in os.listdir(output_dir):
        if not filename.endswith("_description_table.tsv"):
            continue
        path = os.path.join(output_dir, filename)
        df = pd.read_csv(path, sep="\t", dtype=str).fillna("")
        if df.empty or "Accession" not in df.columns:
            continue
        df["Source_File"] = os.path.splitext(filename)[0]
        dataframes.append(df)
    return dataframes


def annotate_accessions(
    config: PipelineConfig,
    all_hits: pd.DataFrame,
    not_found: Set[str],
) -> pd.DataFrame:
    """Fetch and merge detailed annotation data for unique accessions."""
    accessions = sorted({
        acc.split("|")[-1]
        for acc in all_hits["Accession"]
        if isinstance(acc, str) and acc and acc != "NA"
    })
    annotations: List[Dict[str, str]] = []

    for index, accession in enumerate(accessions, start=1):
        print(f"  [{index}/{len(accessions)}] Fetching metadata for {accession}...", end="", flush=True)
        xml_data = entrez_fetch_xml(accession, retries=3)
        if not xml_data:
            print(" FAIL")
            not_found.add(accession)
            annotations.append(
                {
                    "Accession": accession,
                    "Description": "",
                    "Taxonomy": "",
                    "Accession Length": "",
                    "Sequencing Technology": "",
                    "Molecule Type": "",
                    "Topology": "",
                    "Assembly": "",
                    "Release Year": "",
                }
            )
            continue

        try:
            root = ET.fromstring(xml_data)
            gbseq = root.find(".//GBSeq")
            if gbseq is None:
                raise ValueError("No GBSeq in XML")

            description = gbseq.findtext("GBSeq_definition") or ""
            taxonomy = gbseq.findtext("GBSeq_taxonomy") or ""
            accession_length = gbseq.findtext("GBSeq_length") or ""
            molecule_type = gbseq.findtext("GBSeq_moltype") or ""
            topology = gbseq.findtext("GBSeq_topology") or ""
            comment = gbseq.findtext("GBSeq_comment") or ""
            year_val = parse_year_from_xml(gbseq)

            if not taxonomy or year_val is None:
                flat = entrez_fetch_flat(accession, retries=2)
                if flat:
                    if not taxonomy:
                        taxonomy = parse_taxonomy_from_flatfile(flat) or ""
                    if year_val is None:
                        match = re.search(r"Date:\s*(\d{2}-[A-Za-z]{3}-\d{4})", flat)
                        if match:
                            try:
                                year_val = dateparser.parse(match.group(1)).year
                            except Exception:  # pragma: no cover - unexpected date format
                                pass

            technology = detect_sequencing_technology(comment)
            assembly_notes = parse_assembly_lines_from_comment(comment)

            annotations.append(
                {
                    "Accession": accession,
                    "Description": description,
                    "Taxonomy": taxonomy,
                    "Accession Length": accession_length,
                    "Sequencing Technology": technology,
                    "Molecule Type": molecule_type,
                    "Topology": topology,
                    "Assembly": assembly_notes,
                    "Release Year": str(year_val) if year_val is not None else "",
                }
            )
            print(" OK")
        except (ET.ParseError, expat.ExpatError, Exception) as error:
            print(f" PARSE_ERROR: {error}")
            not_found.add(accession)
            annotations.append(
                {
                    "Accession": accession,
                    "Description": "",
                    "Taxonomy": "",
                    "Accession Length": "",
                    "Sequencing Technology": "",
                    "Molecule Type": "",
                    "Topology": "",
                    "Assembly": "",
                    "Release Year": "",
                }
            )
        time.sleep(config.effective_delay)

    annotation_df = pd.DataFrame(annotations)
    all_hits["Accession_clean"] = all_hits["Accession"].astype(str).str.split("|").str[-1]
    annotation_df["Accession_clean"] = annotation_df["Accession"]
    merged = all_hits.merge(
        annotation_df.drop(columns=["Accession"]),
        on="Accession_clean",
        how="left",
        suffixes=("_x", "_y"),
    )
    if "Accession_y" in merged.columns:
        merged.drop(columns=["Accession_y"], inplace=True)
    merged.rename(columns={"Accession_x": "Accession"}, inplace=True)

    for source_file, subset in merged.groupby("Source_File", dropna=False):
        out_tsv = (
            os.path.join(config.output_dir, f"{source_file}.tsv")
            if not source_file.endswith(".tsv")
            else os.path.join(config.output_dir, source_file)
        )
        subset.to_csv(out_tsv, sep="\t", index=False)

    merged.to_pickle(config.partial_pickle)
    return merged


def run_step2(config: PipelineConfig) -> Tuple[pd.DataFrame, Set[str]]:
    """Centralise annotation for unique accessions and persist intermediate state."""
    print(f"\n[STEP 2] Annotating accessions from GenBank — state: {Path(config.partial_pickle).resolve()}")

    not_found: Set[str] = set()
    if os.path.exists(config.partial_pickle) and not config.force_reannotate:
        print(f"[INFO] Resuming from partial state: {config.partial_pickle}")
        all_hits = pd.read_pickle(config.partial_pickle)
        return all_hits, not_found

    dataframes = collect_hit_tables(config.output_dir)
    if not dataframes:
        all_hits = pd.DataFrame(columns=TSV_HEADER + ["Source_File"])
    else:
        all_hits = pd.concat(dataframes, ignore_index=True)

    if all_hits.empty:
        all_hits.to_pickle(config.partial_pickle)
        return all_hits, not_found

    annotated_hits = annotate_accessions(config, all_hits, not_found)
    return annotated_hits, not_found


def run_step3(config: PipelineConfig, all_hits: pd.DataFrame, not_found: Set[str]) -> None:
    """Generate Excel summaries and missing-accession reports."""
    print("\n[STEP 3] Exporting Excel summary")

    with pd.ExcelWriter(config.final_excel, engine="openpyxl") as writer:
        all_hits.to_excel(writer, index=False, sheet_name="Annotations")

        filter_columns = ["Identity(%)", "Coverage(%)", "Align_Length", "Pass_Filter"]
        present_cols = [column for column in filter_columns if column in all_hits.columns]
        if present_cols:
            filtered = all_hits.copy()
            if "Identity(%)" in filtered.columns:
                filtered["Identity(%)"] = pd.to_numeric(filtered["Identity(%)"], errors="coerce")
            if "Coverage(%)" in filtered.columns:
                filtered["Coverage(%)"] = pd.to_numeric(filtered["Coverage(%)"], errors="coerce")
            if "Align_Length" in filtered.columns:
                filtered["Align_Length"] = pd.to_numeric(filtered["Align_Length"], errors="coerce")
            filtered = filtered[
                (filtered.get("Pass_Filter", "") == "YES")
                |
                (
                    (filtered.get("Identity(%)", 0) >= float(config.min_identity))
                    & (filtered.get("Coverage(%)", 0) >= float(config.min_coverage))
                    & (filtered.get("Align_Length", 0) >= int(config.min_matchlen))
                )
            ]
            filtered.to_excel(writer, index=False, sheet_name="Filtered_hits")

        presence = all_hits.groupby(["Accession", "Source_File"]).size().unstack(fill_value=0)
        presence = (presence > 0).replace({True: "✓", False: ""})

        metadata_columns = [
            column
            for column in [
                "Accession",
                "Description",
                "Taxonomy_y",
                "Accession Length",
                "Sequencing Technology",
                "Molecule Type",
                "Topology",
                "Assembly",
                "Release Year_y",
            ]
            if column in all_hits.columns
        ]
        metadata = (
            all_hits[metadata_columns]
            .drop_duplicates(subset="Accession")
            .set_index("Accession")
        )
        matrix_df = metadata.join(presence)
        matrix_df.to_excel(writer, sheet_name="Presence_Matrix")

        if os.path.exists(config.partial_pickle):
            missing_meta = set(
                all_hits.loc[
                    (
                        all_hits.get("Description", "").astype(str) == ""
                        & (all_hits.get("Taxonomy_y", "").astype(str) == "")
                        & (all_hits.get("Accession Length", "").astype(str) == "")
                    ),
                    "Accession",
                ]
                .astype(str)
                .str.split("|")
                .str[-1]
                .unique()
            )
            not_found = not_found.union(missing_meta)

        if not_found:
            pd.DataFrame({"Not_Found": sorted(not_found)}).to_excel(writer, index=False, sheet_name="Not_Found")

    if not_found:
        with open(config.not_found_file, "w", encoding="utf-8") as handle:
            for accession in sorted(not_found):
                handle.write(f"{accession}\n")
        print(f"[INFO] {len(not_found)} accessions failed. Saved to {config.not_found_file}")
    else:
        print("[INFO] All accessions successfully annotated")

    print(f"[DONE] Output saved to: {config.final_excel}")


def run_pipeline(config: PipelineConfig) -> None:
    """Run the full adapter contamination pipeline."""
    configure_logging()
    socket.setdefaulttimeout(30)
    os.makedirs(config.output_dir, exist_ok=True)
    configure_entrez(config)

    minimal_metadata_cache: Dict[str, Dict[str, str]] = {}
    run_step1(config, minimal_metadata_cache)
    all_hits, not_found = run_step2(config)
    run_step3(config, all_hits, not_found)


def build_parser() -> argparse.ArgumentParser:
    """Return the argument parser for the command-line interface."""
    parser = argparse.ArgumentParser(
        description=(
            "Short BLASTn or local search for adapter/barcode sequences, annotate with Entrez, and export Excel."
        )
    )
    parser.add_argument(
        "-i",
        "--input",
        type=str,
        default="ONT_and_Porechop_adapters_and_barcodes_CLEANED.fasta",
        help="Input FASTA with adapter/barcode sequences.",
    )
    parser.add_argument(
        "-o",
        "--output",
        type=str,
        default="shortblast_outputs",
        help="Output directory for per-query TSVs.",
    )
    parser.add_argument(
        "-s",
        "--search_target",
        type=str,
        help="Optional FASTA to search locally instead of BLASTing against GenBank.",
    )
    parser.add_argument(
        "--force-reannotate",
        action="store_true",
        help="Ignore existing partial_annotation.pkl and re-annotate all accessions.",
    )
    parser.add_argument(
        "--email",
        type=str,
        default=os.environ.get("NCBI_EMAIL", ""),
        help="Email for NCBI Entrez (default from env NCBI_EMAIL).",
    )
    parser.add_argument(
        "--api-key",
        type=str,
        default=os.environ.get("NCBI_API_KEY", ""),
        help="NCBI API key (default from env NCBI_API_KEY).",
    )
    parser.add_argument(
        "--min-length",
        type=int,
        default=8,
        help="Minimum query length to process.",
    )
    parser.add_argument(
        "--expect",
        type=float,
        default=10.0,
        help="E-value threshold for BLASTn.",
    )
    parser.add_argument(
        "--hitlist-size",
        type=int,
        default=25,
        help=(
            "Max number of BLAST hits to consider per query (only used by legacy qblast path; ignored in CLI modes)."
        ),
    )
    parser.add_argument(
        "--delay",
        type=float,
        default=None,
        help="Delay between Entrez requests in seconds. Default 0.34 if API key present, else 0.5.",
    )
    parser.add_argument(
        "--min-identity",
        type=float,
        default=95.0,
        help="Minimum percent identity to keep a hit (default 95).",
    )
    parser.add_argument(
        "--min-coverage",
        type=float,
        default=80.0,
        help="Minimum percent coverage of the QUERY to keep a hit (default 80).",
    )
    parser.add_argument(
        "--min-matchlen",
        type=int,
        default=16,
        help="Minimum aligned length in bp to keep a hit (default 16).",
    )
    parser.add_argument(
        "--blast-mode",
        type=str,
        default="local",
        choices=["local", "remote"],
        help="Where to run BLAST: 'local' uses local blastn against --db; 'remote' uses blastn -remote against NCBI. Default: local.",
    )
    parser.add_argument(
        "--blastn",
        type=str,
        default="blastn",
        help="Path to blastn executable (default: blastn in PATH).",
    )
    parser.add_argument(
        "--db",
        type=str,
        default="nt",
        help="BLAST database to use (local DB name/path for --blast-mode local; 'nt' for remote).",
    )
    parser.add_argument(
        "--max-target-seqs",
        type=int,
        default=100_000_000,
        help="Very high cap to avoid truncation. Used by blastn CLI.",
    )
    return parser


__all__ = [
    "PipelineConfig",
    "run_pipeline",
    "build_parser",
]
