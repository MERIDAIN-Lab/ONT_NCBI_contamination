#!/usr/bin/env python3
# Combined BLAST + Entrez annotation with robust fallbacks
# - STEP 1: For each query in FASTA, run BLAST (or local search) and write a TSV per query.
#           Each TSV row includes Organism, Taxonomy, and Release Year fetched via Entrez (with flatfile fallback).
# - STEP 2: Centralized re-annotation of unique accessions (definition, taxonomy, accession length, molecule type,
#           topology, sequencing technology, assembly notes, release year). Merge back to produce _x/_y columns.
# - STEP 3: Export Excel with Annotations, Presence_Matrix, Not_Found. Idempotent; supports resume and force-reannotate.

import os
import re
import time
import socket
import logging
import argparse
import subprocess, tempfile
from io import StringIO
from contextlib import closing
from typing import Dict, List, Optional, Tuple

import pandas as pd
from dateutil import parser as dateparser
from Bio import SeqIO, Entrez
import xml.etree.ElementTree as ET
import xml.parsers.expat as expat

# Optional: pairwise local alignment for --search_target mode
try:
    from Bio import pairwise2
    HAVE_PAIRWISE2 = True
except Exception:
    HAVE_PAIRWISE2 = False

# ---------------------------------------------------------
# Args
# ---------------------------------------------------------
parser = argparse.ArgumentParser(
    description="Short BLASTn or local search for adapter/barcode sequences, annotate with Entrez, and export Excel."
)
parser.add_argument(
    "-i", "--input", type=str, default="ONT_and_Porechop_adapters_and_barcodes_CLEANED.fasta",
    help="Input FASTA with adapter/barcode sequences."
)
parser.add_argument(
    "-o", "--output", type=str, default="shortblast_outputs",
    help="Output directory for per-query TSVs."
)
parser.add_argument(
    "-s", "--search_target", type=str,
    help="Optional FASTA to search locally instead of BLASTing against GenBank."
)
parser.add_argument(
    "--force-reannotate", action="store_true",
    help="Ignore existing partial_annotation.pkl and re-annotate all accessions."
)
parser.add_argument(
    "--email", type=str, default=os.environ.get("NCBI_EMAIL", ""),
    help="Email for NCBI Entrez (default from env NCBI_EMAIL)."
)
parser.add_argument(
    "--api-key", type=str, default=os.environ.get("NCBI_API_KEY", ""),
    help="NCBI API key (default from env NCBI_API_KEY)."
)
parser.add_argument(
    "--min-length", type=int, default=8,
    help="Minimum query length to process."
)
parser.add_argument(
    "--expect", type=float, default=10.0,
    help="E-value threshold for BLASTn."
)
parser.add_argument(
    "--hitlist-size", type=int, default=25,
    help="Max number of BLAST hits to consider per query (only used by legacy qblast path; ignored in CLI modes)."
)
parser.add_argument(
    "--delay", type=float, default=None,
    help="Delay between Entrez requests in seconds. Default 0.34 if API key present, else 0.5."
)
parser.add_argument("--min-identity", type=float, default=95.0,
    help="Minimum percent identity to keep a hit (default 95).")
parser.add_argument("--min-coverage", type=float, default=80.0,
    help="Minimum percent coverage of the QUERY to keep a hit (default 80).")
parser.add_argument("--min-matchlen", type=int, default=16,
    help="Minimum aligned length in bp to keep a hit (default 16).")

parser.add_argument(
    "--blast-mode", type=str, default="local", choices=["local", "remote"],
    help="Where to run BLAST: 'local' uses local blastn against --db; 'remote' uses blastn -remote against NCBI. Default: local."
)
parser.add_argument(
    "--blastn", type=str, default="blastn",
    help="Path to blastn executable (default: blastn in PATH)."
)
parser.add_argument(
    "--db", type=str, default="nt",
    help="BLAST database to use (local DB name/path for --blast-mode local; 'nt' for remote)."
)
parser.add_argument(
    "--max-target-seqs", type=int, default=100000000,
    help="Very high cap to avoid truncation. Used by blastn CLI."
)




args = parser.parse_args()

# ---------------------------------------------------------
# Config
# ---------------------------------------------------------
socket.setdefaulttimeout(30)
logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")

OUTPUT_DIR = args.output
INPUT_FASTA = args.input
SEARCH_TARGET = args.search_target
MIN_LENGTH = int(args.min_length)
FINAL_EXCEL = "Combined_BLAST_Report.xlsx"
NOT_FOUND_FILE = "accessions_not_retrieved.txt"
PARTIAL_PKL = "partial_annotation.pkl"

os.makedirs(OUTPUT_DIR, exist_ok=True)

# Entrez
Entrez.email = args.email or "your.email@example.com"
if args.api_key:
    Entrez.api_key = args.api_key
DELAY = args.delay if args.delay is not None else (0.34 if args.api_key else 0.5)

# ---------------------------------------------------------
# Helpers
# ---------------------------------------------------------

def run_blast_cli_return_rows(seq: str) -> List[List[str]]:
    """
    Run blastn (local or remote) for a short query and return rows matching TSV_HEADER order.
    Uses -outfmt 6 with the following columns:
    sacc evalue pident length qlen qstart qend sstart send qcovs slen
    """
    # Write query to a temp FASTA
    with tempfile.NamedTemporaryFile(mode="w", suffix=".fa", delete=False) as tf:
        tf.write(">q\n")
        tf.write(seq + "\n")
        qpath = tf.name
    cmd = [
        args.blastn,
        "-task", "blastn-short",
        "-query", qpath,
        "-word_size", "7",
        "-evalue", str(args.expect),
        "-dust", "no",
        "-soft_masking", "false",
        "-db", args.db,
        "-outfmt", "6 sacc evalue pident length qlen qstart qend sstart send qcovs slen",
        "-max_target_seqs", str(args.max_target_seqs),
    ]
    if args.blast_mode == "remote":
        cmd.append("-remote")
    # Execute
    try:
        res = subprocess.run(cmd, check=True, capture_output=True, text=True)
        lines = [ln for ln in res.stdout.splitlines() if ln.strip()]
    except subprocess.CalledProcessError as e:
        logging.error(f"blastn failed: {e}\nSTDERR: {e.stderr}")
        lines = []
    finally:
        try:
            os.unlink(qpath)
        except Exception:
            pass

    rows: List[List[str]] = []
    for ln in lines:
        parts = ln.split("\t")
        if len(parts) < 11:
            continue
        (sacc, evalue, pident, alen, qlen, qstart, qend, sstart, send, qcovs, slen) = parts[:11]
        # Minimal meta now; full meta added later in STEP 1 via get_minimal_meta
        meta = get_minimal_meta(sacc)
        organism = meta.get("Organism") or "NA"
        taxonomy = meta.get("Taxonomy", "")
        year = meta.get("Release Year", "")
        # Filters
        try:
            id_pct = float(pident)
            cov_pct = float(qcovs)
            alen_i = int(float(alen))
        except Exception:
            id_pct, cov_pct, alen_i = 0.0, 0.0, 0
        pass_filter = (
            (alen_i >= int(args.min_matchlen)) and
            (id_pct >= float(args.min_identity)) and
            (cov_pct >= float(args.min_coverage))
        )
        row = [
            sacc,            # Accession
            slen,            # Length (subject length)
            evalue,          # E-value
            "",              # Definition (fill later if desired)
            organism,        # Organism
            taxonomy,        # Taxonomy
            year,            # Release Year
            qstart, qend, sstart, send,  # coords
            alen,            # Align_Length
            qlen,            # Query_Length
            f"{id_pct:.2f}", # Identity(%)
            f"{cov_pct:.2f}",# Coverage(%)
            "YES" if pass_filter else "NO"
        ]
        rows.append(row)
    return rows

TSV_HEADER = [
    "Accession", "Length", "E-value", "Definition",
    "Organism", "Taxonomy", "Release Year",
    "query_start", "query_end", "subject_start", "subject_end",
    "Align_Length", "Query_Length", "Identity(%)", "Coverage(%)", "Pass_Filter"
]

def sanitize_id(s: str, maxlen: int = 64) -> str:
    clean = re.sub(r"[^A-Za-z0-9]+", "_", s).strip("_")
    return clean[:maxlen] if len(clean) > maxlen else clean

def tsv_path_for(seq_id: str, seq: str) -> str:
    clean_id = sanitize_id(seq_id)
    suffix_core = seq if len(seq) <= 30 else f"{seq[:15]}...{seq[-10:]}"
    suffix = sanitize_id(suffix_core, maxlen=40)
    fname = f"{clean_id}__{suffix}_description_table.tsv"
    return os.path.join(OUTPUT_DIR, fname)

def write_tsv_header(path: str) -> None:
    with open(path, "w") as out:
        out.write("\t".join(TSV_HEADER) + "\n")

def append_tsv_row(path: str, row: List[str]) -> None:
    with open(path, "a") as out:
        out.write("\t".join(row) + "\n")

def select_best_hsp(alignment) -> Optional[Tuple]:
    # Best by e-value, then bits (desc), then align length (desc)
    best = None
    for hsp in alignment.hsps:
        e = getattr(hsp, "expect", float("inf"))
        bits = getattr(hsp, "bits", 0.0)
        alen = getattr(hsp, "align_length", getattr(hsp, "align_len", 0))
        key = (e, -bits, -alen)
        if best is None or key < best[1]:
            best = (hsp, key)
    return best[0] if best else None

def parse_accession(hit_id: str) -> str:
    m = re.search(r"\|(ref|gb|dbj|emb)\|([^\|]+)\|", hit_id)
    if m:
        return m.group(2)
    tokens = hit_id.split()
    return tokens[0]

def entrez_fetch_xml(acc: str, retries: int = 3) -> Optional[str]:
    for attempt in range(1, retries + 1):
        try:
            with closing(Entrez.efetch(db="nuccore", id=acc, retmode="xml")) as handle:
                return handle.read()
        except Exception as e:
            logging.warning(f"efetch XML failed for {acc} (attempt {attempt}): {e}")
            time.sleep(1.5 * attempt)
    return None

def entrez_fetch_flat(acc: str, retries: int = 2) -> Optional[str]:
    for attempt in range(1, retries + 1):
        try:
            with closing(Entrez.efetch(db="nuccore", id=acc, rettype="gb", retmode="text")) as handle:
                return handle.read()
        except Exception as e:
            logging.warning(f"efetch flatfile failed for {acc} (attempt {attempt}): {e}")
            time.sleep(1.5 * attempt)
    return None

def parse_year_from_xml(root_el) -> Optional[int]:
    dt = root_el.findtext(".//GBSeq_update-date") or root_el.findtext(".//GBSeq_create-date")
    if not dt:
        return None
    try:
        return dateparser.parse(dt).year
    except Exception:
        m = re.search(r"\d{2}-[A-Za-z]{3}-(\d{4})", dt)
        if m:
            try:
                return int(m.group(1))
            except Exception:
                return None
    return None

def parse_taxonomy_from_flatfile(gb_text: str) -> str:
    # Extract the multiline taxonomy block that follows the ORGANISM line
    if not gb_text:
        return ""
    m = re.search(r"^ORGANISM[ \t]+[^\n]*\n((?:[ \t]{2,}.+\n)+)", gb_text, flags=re.MULTILINE)
    if not m:
        # Some records include 'TAXONOMY' lines in other formats
        m2 = re.search(r"^TAXONOMY[ \t]+(.+)$", gb_text, flags=re.MULTILINE)
        if m2:
            return m2.group(1).strip().rstrip(".")
        return ""
    block = m.group(1)
    lines = [ln.strip() for ln in block.splitlines()]
    lineage = " ".join(lines)
    # GenBank lineage already semicolon-delimited; remove trailing period if present
    lineage = lineage.rstrip(".")
    return lineage

def parse_organism_from_flatfile(gb_text: str) -> str:
    if not gb_text:
        return "NA"
    m = re.search(r"^ORGANISM[ \t]+(.+)$", gb_text, flags=re.MULTILINE)
    return m.group(1).strip() if m else "NA"

def parse_assembly_lines_from_comment(comment_text: str) -> str:
    if not comment_text:
        return ""
    lines = [ln.strip() for ln in comment_text.splitlines() if "assembly" in ln.lower()]
    return "; ".join(lines)

def detect_sequencing_technology(comment_text: str) -> str:
    if not comment_text:
        return ""
    c = comment_text.lower()
    nanopore_keywords = ["nanopore", "ont", "minion", "gridion", "promethion", "flongle", "minknow", "oxford nanopore"]
    if any(k in c for k in nanopore_keywords):
        return "Oxford Nanopore"
    return ""

# Minimal per-accession cache used in STEP 1 to avoid repeated calls
min_meta_cache: Dict[str, Dict[str, str]] = {}

def get_minimal_meta(acc: str) -> Dict[str, str]:
    if acc in min_meta_cache:
        return min_meta_cache[acc]
    meta = {"Organism": "NA", "Taxonomy": "", "Release Year": ""}
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
        except (ET.ParseError, expat.ExpatError) as e:
            logging.warning(f"XML parse error for {acc}: {e}")

    # Fallbacks if missing
    if not meta["Taxonomy"] or meta["Organism"] == "NA" or not meta["Release Year"]:
        flat = entrez_fetch_flat(acc, retries=2)
        if flat:
            if not meta["Organism"] or meta["Organism"] == "NA":
                meta["Organism"] = parse_organism_from_flatfile(flat) or "NA"
            if not meta["Taxonomy"]:
                meta["Taxonomy"] = parse_taxonomy_from_flatfile(flat) or ""
            if not meta["Release Year"]:
                ym = re.search(r"Date:\s*(\d{2}-[A-Za-z]{3}-\d{4})", flat)
                if ym:
                    try:
                        meta["Release Year"] = str(dateparser.parse(ym.group(1)).year)
                    except Exception:
                        pass
    min_meta_cache[acc] = meta
    time.sleep(DELAY)
    return meta

# ---------------------------------------------------------
# STEP 1: Generate TSVs via BLAST or local search and include Organism/Taxonomy/Year
# ---------------------------------------------------------
print("[STEP 1] Generating TSVs per query (skip existing)")

if SEARCH_TARGET:
    if not HAVE_PAIRWISE2:
        raise RuntimeError("pairwise2 not available. Install Biopython to use --search_target.")
    target_seqs = {rec.id: str(rec.seq) for rec in SeqIO.parse(SEARCH_TARGET, "fasta")}
else:
    target_seqs = {}

for rec in SeqIO.parse(INPUT_FASTA, "fasta"):
    seq_id = rec.id
    seq = str(rec.seq)
    if len(seq) < MIN_LENGTH:
        print(f"[SKIP] {seq_id}: sequence too short ({len(seq)} bp)")
        continue

    out_path = tsv_path_for(seq_id, seq)
    if os.path.exists(out_path):
        print(f"[SKIP] {seq_id}: already processed -> {os.path.basename(out_path)}")
        continue

    print(f"[INFO] Processing {seq_id} ({len(seq)} bp) -> {os.path.basename(out_path)}")
    write_tsv_header(out_path)

    if SEARCH_TARGET:
        # Local best-match search across target FASTA
        best = None  # (target_id, score)
        for tid, tseq in target_seqs.items():
            try:
                alns = pairwise2.align.localms(seq, tseq, 2, -1, -0.5, -0.1)
            except Exception as e:
                logging.warning(f"Local align failed for {seq_id} vs {tid}: {e}")
                alns = []
            if alns:
                score = alns[0].score
                if best is None or score > best[1]:
                    best = (tid, score)
        if best:
            # Local mode cannot Entrez-annotate arbitrary target IDs reliably; leave blanks for Taxonomy/Year
            row = [
                best[0], str(len(target_seqs[best[0]])), "NA", "Local_match",
                "NA", "", "", "NA", "NA", "NA", "NA",
                "NA", str(len(seq)), "NA", "NA", "NA"
            ]
            append_tsv_row(out_path, row)
        else:
            append_tsv_row(out_path, ["NA","NA","NA","No_match_found","NA","","","NA","NA","NA","NA","NA",str(len(seq)),"NA","NA","NA"])
        print(f"[OK] Local search saved: {os.path.basename(out_path)}")
        continue

    # BLASTn short via CLI (local or remote) with effectively no truncation
    rows = run_blast_cli_return_rows(seq)
    if not rows:
        append_tsv_row(out_path, [
            "NA","NA","NA","No_hits","NA","","","NA","NA","NA","NA","NA",
            str(len(seq)), "NA","NA","NA"
        ])
    else:
        for row in rows:
            append_tsv_row(out_path, row)
    print(f"[OK] BLAST results saved: {os.path.basename(out_path)}")

# ---------------------------------------------------------
# STEP 2: Centralized annotation of unique accessions and merge (_x/_y columns)
# ---------------------------------------------------------
print(f"\n[STEP 2] Annotating accessions from GenBank — state: {os.path.abspath(PARTIAL_PKL)}")

not_found: set = set()
if os.path.exists(PARTIAL_PKL) and not args.force_reannotate:
    print(f"[INFO] Resuming from partial state: {PARTIAL_PKL}")
    all_hits = pd.read_pickle(PARTIAL_PKL)
else:
    dfs: List[pd.DataFrame] = []
    for fname in os.listdir(OUTPUT_DIR):
        if not fname.endswith("_description_table.tsv"):
            continue
        path = os.path.join(OUTPUT_DIR, fname)
        df = pd.read_csv(path, sep="\t", dtype=str).fillna("")
        if df.empty or "Accession" not in df.columns:
            continue
        df["Source_File"] = os.path.splitext(fname)[0]
        dfs.append(df)
    if not dfs:
        all_hits = pd.DataFrame(columns=TSV_HEADER + ["Source_File"])
    else:
        all_hits = pd.concat(dfs, ignore_index=True)

    # Unique valid accessions
    accs = sorted({a.split("|")[-1] for a in all_hits["Accession"] if isinstance(a, str) and a and a != "NA"})
    annotations: List[Dict[str, str]] = []

    for i, acc in enumerate(accs, 1):
        print(f"  [{i}/{len(accs)}] Fetching metadata for {acc}...", end="", flush=True)
        xml_data = entrez_fetch_xml(acc, retries=3)
        if not xml_data:
            print(" FAIL")
            not_found.add(acc)
            # Fill with blanks to keep schema
            annotations.append({
                "Accession": acc,
                "Description": "",
                "Taxonomy": "",
                "Accession Length": "",
                "Sequencing Technology": "",
                "Molecule Type": "",
                "Topology": "",
                "Assembly": "",
                "Release Year": ""
            })
            continue

        try:
            root = ET.fromstring(xml_data)
            gbseq = root.find(".//GBSeq")
            if gbseq is None:
                raise ValueError("No GBSeq in XML")

            description = gbseq.findtext("GBSeq_definition") or ""
            taxonomy   = gbseq.findtext("GBSeq_taxonomy") or ""
            acc_len    = gbseq.findtext("GBSeq_length") or ""
            moltype    = gbseq.findtext("GBSeq_moltype") or ""
            topology   = gbseq.findtext("GBSeq_topology") or ""
            comment    = gbseq.findtext("GBSeq_comment") or ""
            year_val   = parse_year_from_xml(gbseq)

            # Fallbacks if missing
            if not taxonomy or year_val is None:
                flat = entrez_fetch_flat(acc, retries=2)
                if flat:
                    if not taxonomy:
                        taxonomy = parse_taxonomy_from_flatfile(flat) or ""
                    if year_val is None:
                        ym = re.search(r"Date:\s*(\d{2}-[A-Za-z]{3}-\d{4})", flat)
                        if ym:
                            try:
                                year_val = dateparser.parse(ym.group(1)).year
                            except Exception:
                                pass

            tech = detect_sequencing_technology(comment)
            assembly_notes = parse_assembly_lines_from_comment(comment)

            annotations.append({
                "Accession": acc,
                "Description": description,
                "Taxonomy": taxonomy,
                "Accession Length": acc_len,
                "Sequencing Technology": tech,
                "Molecule Type": moltype,
                "Topology": topology,
                "Assembly": assembly_notes,
                "Release Year": str(year_val) if year_val is not None else ""
            })
            print(" OK")
        except (ET.ParseError, expat.ExpatError, Exception) as e:
            print(f" PARSE_ERROR: {e}")
            not_found.add(acc)
            annotations.append({
                "Accession": acc,
                "Description": "",
                "Taxonomy": "",
                "Accession Length": "",
                "Sequencing Technology": "",
                "Molecule Type": "",
                "Topology": "",
                "Assembly": "",
                "Release Year": ""
            })
        time.sleep(DELAY)

    ann_df = pd.DataFrame(annotations)

    # Merge to produce duplicated columns with suffixes (_x from STEP 1, _y from STEP 2)
    # Join on cleaned accession to be safe
    all_hits["Accession_clean"] = all_hits["Accession"].astype(str).str.split("|").str[-1]
    ann_df["Accession_clean"] = ann_df["Accession"]
    all_hits = all_hits.merge(
        ann_df.drop(columns=["Accession"]),
        on="Accession_clean",
        how="left",
        suffixes=("_x", "_y")
    )
    # Keep original Accession from STEP 1
    if "Accession_y" in all_hits.columns:
        all_hits.drop(columns=["Accession_y"], inplace=True)
    all_hits.rename(columns={"Accession_x": "Accession"}, inplace=True)

    # Persist enriched TSVs back to disk for each Source_File
    for sf, sub in all_hits.groupby("Source_File", dropna=False):
        out_tsv = os.path.join(OUTPUT_DIR, f"{sf}.tsv") if not sf.endswith(".tsv") else os.path.join(OUTPUT_DIR, sf)
        sub.to_csv(out_tsv, sep="\t", index=False)

    # Persist pickle for quick resume
    all_hits.to_pickle(PARTIAL_PKL)

# ---------------------------------------------------------
# STEP 3: Excel export
# ---------------------------------------------------------
print("\n[STEP 3] Exporting Excel summary")

with pd.ExcelWriter(FINAL_EXCEL, engine="openpyxl") as writer:
    # Full annotations with _x/_y columns preserved
    all_hits.to_excel(writer, index=False, sheet_name="Annotations")



    flt_cols = ["Identity(%)", "Coverage(%)", "Align_Length", "Pass_Filter"]
    present_cols = [c for c in flt_cols if c in all_hits.columns]
    if present_cols:
        # cast to float de forma segura
        ah = all_hits.copy()
        for c in ["Identity(%)", "Coverage(%)"]:
            if c in ah.columns:
                ah[c] = pd.to_numeric(ah[c], errors="coerce")
        if "Align_Length" in ah.columns:
            ah["Align_Length"] = pd.to_numeric(ah["Align_Length"], errors="coerce")
        filtered = ah[
            (ah.get("Pass_Filter", "") == "YES") |
            (
                (ah.get("Identity(%)", 0) >= float(args.min_identity)) &
                (ah.get("Coverage(%)", 0) >= float(args.min_coverage)) &
                (ah.get("Align_Length", 0) >= int(args.min_matchlen))
            )
        ]
        filtered.to_excel(writer, index=False, sheet_name="Filtered_hits")




    

    # Presence matrix
    presence = all_hits.groupby(["Accession", "Source_File"]).size().unstack(fill_value=0)
    presence = (presence > 0).replace({True: "✓", False: ""})

    # Metadata view keyed by Accession (take STEP 2 fields)
    meta_cols = [c for c in [
        "Accession", "Description", "Taxonomy_y", "Accession Length",
        "Sequencing Technology", "Molecule Type", "Topology", "Assembly", "Release Year_y"
    ] if c in all_hits.columns]
    metadata = all_hits[meta_cols].drop_duplicates(subset="Accession").set_index("Accession")
    matrix_df = metadata.join(presence)
    matrix_df.to_excel(writer, sheet_name="Presence_Matrix")

    # Not Found sheet
    if os.path.exists(PARTIAL_PKL):
        missing_meta = set(
            all_hits.loc[
                (all_hits.get("Description", "").astype(str) == "") &
                (all_hits.get("Taxonomy_y", "").astype(str) == "") &
                (all_hits.get("Accession Length", "").astype(str) == ""),
                "Accession"
            ].astype(str).str.split("|").str[-1].unique()
        )
        not_found = not_found.union(missing_meta)

    if not_found:
        pd.DataFrame({"Not_Found": sorted(not_found)}).to_excel(writer, index=False, sheet_name="Not_Found")

# Also write a plain text list
if not_found:
    with open(NOT_FOUND_FILE, "w") as f:
        for acc in sorted(not_found):
            f.write(acc + "\n")
    print(f"[INFO] {len(not_found)} accessions failed. Saved to", NOT_FOUND_FILE)
else:
    print("[INFO] All accessions successfully annotated")

print(f"[DONE] Output saved to: {FINAL_EXCEL}")
