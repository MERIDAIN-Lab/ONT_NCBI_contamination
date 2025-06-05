import os
import re
import time
import logging
import pandas as pd
from Bio import SeqIO, Entrez
from Bio.Blast import NCBIWWW, NCBIXML
from io import StringIO
from openpyxl import load_workbook
import xml.etree.ElementTree as ET
import urllib.error
import argparse

# ─── ARGUMENT PARSING ─────────────────────────────────────────────
parser = argparse.ArgumentParser(
    description="Run BLASTn annotation on a FASTA file or use GenBank accessions to annotate previously generated hits."
)
parser.add_argument(
    "-i", "--input", type=str, default="ONT_and_Porechop_adapters_and_barcodes_CLEANED.fasta",
    help="Input FASTA file with adapter/barcode sequences (default: ONT_and_Porechop_adapters_and_barcodes_CLEANED.fasta)"
)
parser.add_argument(
    "-o", "--output", type=str, default="shortblast_outputs",
    help="Output directory where BLAST and annotation results will be stored (default: shortblast_outputs)"
)
parser.add_argument(
    "-s", "--search_target", type=str,
    help="Optional FASTA file to search for adapters locally instead of BLASTing against GenBank"
)
args = parser.parse_args()

# Assign search_target argument
search_target = args.search_target

# ─── CONFIGURATION ────────────────────────────────────────────────
Entrez.email = "florenmartino@gmail.com"
Entrez.api_key = "36bf13267a495a3f7007c3a8d386587e0308"

input_fasta = args.input
output_dir = args.output
min_length = 8
final_excel = "Combined_BLAST_Report.xlsx"
matrix_file = "Accession_vs_Adapters_Matrix.xlsx"
not_found_file = "accessions_not_retrieved.txt"
partial_pickle = "partial_annotation.pkl"
delay_between_requests = 0.34

os.makedirs(output_dir, exist_ok=True)

# ─── PART 1: BLASTn SHORT & SAVE INDIVIDUAL TSV FILES ─────────────
print("[STEP 1] Running BLASTn and saving .tsv tables (skip existing)")

for record in SeqIO.parse(input_fasta, "fasta"):
    seq_id = record.id
    seq = str(record.seq)
    clean_id = re.sub(r"[^a-zA-Z0-9]+", "_", seq_id)
    suffix = f"_{seq}" if len(seq) < 100 else f"_{seq[:50]}[...]"
    table_filename = f"{clean_id}{suffix}_description_table.tsv"
    table_path = os.path.join(output_dir, table_filename)

    if os.path.exists(table_path):
        print(f"[SKIP] {seq_id}: already processed")
        continue

    if len(seq) < min_length:
        print(f"[SKIP] {seq_id}: sequence too short ({len(seq)} bp)")
        continue

    print(f"[INFO] BLASTing {seq_id} ({len(seq)} bp)...")
    try:
        if search_target:
            from Bio import pairwise2
            target_seqs = {rec.id: str(rec.seq) for rec in SeqIO.parse(search_target, "fasta")}
            best_hit = None
            best_score = 0
            for target_id, target_seq in target_seqs.items():
                alignments = pairwise2.align.localms(seq, target_seq, 2, -1, -0.5, -0.1)
                if alignments:
                    score = alignments[0].score
                    if score > best_score:
                        best_score = score
                        best_hit = (target_id, score, alignments[0])
            with open(table_path, "w") as out:
                out.write("Accession\tLength\tE-value\tDefinition\tOrganism\tquery_start\tquery_end\tsubject_start\tsubject_end\n")
                if best_hit:
                    aln = best_hit[2]
                    out.write(f"{best_hit[0]}\t{len(aln.seqB)}\tNA\tLocal_match\tNA\t{aln.start}\t{aln.end}\tNA\tNA\n")
                else:
                    out.write("NA\tNA\tNA\tNo_match_found\tNA\tNA\tNA\tNA\tNA\n")
            print(f"[OK] Local search saved: {table_path}")
        else:
            result_handle = NCBIWWW.qblast(
                program="blastn",
                database="nt",
                sequence=seq,
                expect=1000,
                word_size=7,
                short_query=True,
                hitlist_size=50,
                format_type="XML"
            )
            blast_data = result_handle.read()
            result_handle.close()
            blast_record = NCBIXML.read(StringIO(blast_data))

            with open(table_path, "w") as out:
                out.write("Accession\tLength\tE-value\tDefinition\tOrganism\tquery_start\tquery_end\tsubject_start\tsubject_end\n")
                for alignment in blast_record.alignments:
                    for hsp in alignment.hsps:
                        acc_match = re.search(r"\|(ref|gb|dbj|emb)\|([^\|]+)\|", alignment.hit_id)
                        accession = acc_match.group(2) if acc_match else alignment.hit_id
                        title = alignment.title
                        definition, organism = (title.rsplit(" [", 1) + ["NA"])[:2]
                        organism = organism.rstrip("]")

                        out.write(f"{accession}\t{alignment.length}\t{hsp.expect:.2e}\t{definition}\t{organism}\t{hsp.query_start}\t{hsp.query_end}\t{hsp.sbjct_start}\t{hsp.sbjct_end}\n")
                        break
            print(f"[OK] Saved: {table_path}")
    except Exception as e:
        print(f"[ERROR] {seq_id}: {e}")

# ─── PART 2: GENBANK METADATA RETRIEVAL ───────────────────────────
print(f"\n[STEP 2] Annotating accessions from GenBank — saving state to: {os.path.abspath(partial_pickle)}")
not_found = set()
all_hits = []

accession_cache = {}

if os.path.exists(partial_pickle):
    print(f"[INFO] Resuming from previous partial state: {partial_pickle}")
    all_hits = pd.read_pickle(partial_pickle)
else:
    for file in os.listdir(output_dir):
        if not file.endswith("_description_table.tsv"):
            continue

        path = os.path.join(output_dir, file)
        df = pd.read_csv(path, sep="\t")
        if df.empty or "Accession" not in df.columns:
            continue

        df["Source_File"] = os.path.splitext(file)[0]
        print(f"[INFO] Annotating {file} ({len(df)} hits)")

        annotations = []

        for i, acc in enumerate(df["Accession"], 1):
            acc_clean = acc.split("|")[-1]
            if acc_clean in accession_cache:
                annotations.append(accession_cache[acc_clean])
                print(f"  [{i}/{len(df)}] Using cached metadata for {acc_clean}")
                continue
            print(f"  [{i}/{len(df)}] Fetching metadata for {acc_clean}...", end="", flush=True)

            max_retries = 3
            for attempt in range(max_retries):
                try:
                    handle = Entrez.efetch(db="nucleotide", id=acc_clean, retmode="xml")
                    xml_data = handle.read()
                    handle.close()
                    break  # successful fetch
                except urllib.error.HTTPError as e:
                    print(f" RETRY {attempt+1} (HTTPError: {e.code})", end="", flush=True)
                    time.sleep(2 * (attempt + 1))
                except Exception as e:
                    print(f" RETRY {attempt+1} (Error: {e})", end="", flush=True)
                    time.sleep(2 * (attempt + 1))
            else:
                print(" FETCH ERROR")
                not_found.add(acc_clean)
                annotations.append({
                    "Accession": acc,
                    "Description": "",
                    "Taxonomy": "",
                    "Accession Length": "",
                    "Sequencing Technology": "",
                    "Molecule Type": "",
                    "Topology": "",
                    "Assembly": ""
                })
                continue

            try:
                root = ET.fromstring(xml_data)
                gbseq = root.find(".//GBSeq")
                if gbseq is None:
                    raise ValueError("No GBSeq found in XML response")

                def get_text(tag):
                    el = gbseq.find(tag)
                    return el.text if el is not None else ""

                comment = get_text("GBSeq_comment").lower()
                tech = "MinION" if any(k in comment for k in ["nanopore", "ont", "minion"]) else ""
                if "sequencing technology:" in comment:
                    tech = comment.split("sequencing technology:")[1].split(";")[0].strip()

                metadata = {
                    "Accession": acc,
                    "Description": get_text("GBSeq_definition"),
                    "Taxonomy": get_text("GBSeq_taxonomy"),
                    "Accession Length": get_text("GBSeq_length"),
                    "Sequencing Technology": tech,
                    "Molecule Type": get_text("GBSeq_moltype"),
                    "Topology": get_text("GBSeq_topology"),
                    "Assembly": get_text("GBSeq_comment") if "assembly" in comment else ""
                }
                accession_cache[acc_clean] = metadata
                annotations.append(metadata)
                print(" OK")

            except Exception as parse_err:
                print(f" PARSE ERROR: {parse_err}")
                not_found.add(acc_clean)
                annotations.append({
                    "Accession": acc,
                    "Description": "",
                    "Taxonomy": "",
                    "Accession Length": "",
                    "Sequencing Technology": "",
                    "Molecule Type": "",
                    "Topology": "",
                    "Assembly": ""
                })

            time.sleep(delay_between_requests)

        ann_df = pd.DataFrame(annotations)
        df = df.merge(ann_df, on="Accession", how="left")
        all_hits.append(df)

    all_hits = pd.concat(all_hits, ignore_index=True)
    all_hits.to_pickle(partial_pickle)

# ─── PART 3: COMBINE ALL AND EXPORT ──────────────────────────────
print("\n[STEP 3] Exporting Excel summary")

with pd.ExcelWriter(final_excel, engine="openpyxl") as writer:
    all_hits.to_excel(writer, index=False, sheet_name="Annotations")

    presence = all_hits.groupby(["Accession", "Source_File"]).size().unstack(fill_value=0)
    presence = presence.applymap(lambda x: "✓" if x > 0 else "")

    meta_cols = ["Accession", "Description", "Taxonomy", "Accession Length",
                 "Sequencing Technology", "Molecule Type", "Topology", "Assembly"]
    metadata = all_hits[meta_cols].drop_duplicates(subset="Accession").set_index("Accession")
    matrix_df = metadata.join(presence)
    matrix_df.to_excel(writer, sheet_name="Presence_Matrix")

    if not_found:
        pd.DataFrame({"Not_Found": sorted(not_found)}).to_excel(writer, index=False, sheet_name="Not_Found")

if not_found:
    with open(not_found_file, "w") as f:
        for acc in sorted(not_found):
            f.write(acc + "\n")
    print(f"[INFO] {len(not_found)} accessions failed. Saved to {not_found_file}")
else:
    print("[INFO] All accessions successfully annotated")

print(f"[DONE] Output saved to: {final_excel}")
