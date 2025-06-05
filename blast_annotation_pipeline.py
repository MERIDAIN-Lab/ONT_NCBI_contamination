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

# ─── CONFIGURATION ────────────────────────────────────────────────
# In local file it was the version "1.2_New_search_complete.py"
Entrez.email = "florenmartino@gmail.com"
Entrez.api_key = "36bf13267a495a3f7007c3a8d386587e0308"

input_fasta = "ONT_and_Porechop_adapters_and_barcodes_CLEANED.fasta"
output_dir = "shortblast_outputs"
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
    table_path = os.path.join(output_dir, f"{seq_id}_description_table.tsv")

    if os.path.exists(table_path):
        print(f"[SKIP] {seq_id}: already processed")
        continue

    if len(seq) < min_length:
        print(f"[SKIP] {seq_id}: sequence too short ({len(seq)} bp)")
        continue

    print(f"[INFO] BLASTing {seq_id} ({len(seq)} bp)...")
    try:
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
print("\n[STEP 2] Annotating accessions from GenBank (skip existing)")
not_found = set()
all_hits = []

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

        df["Source_File"] = file.replace("_description_table.tsv", "")
        print(f"[INFO] Annotating {file} ({len(df)} hits)")

        annotations = []

        for i, acc in enumerate(df["Accession"], 1):
            acc_clean = acc.split("|")[-1]
            print(f"  [{i}/{len(df)}] Fetching metadata for {acc_clean}...", end="")
            try:
                handle = Entrez.efetch(db="nuccore", id=acc_clean, retmode="xml")
                xml_data = handle.read()
                handle.close()
                root = ET.fromstring(xml_data)
                gbseq = root.find(".//GBSeq")

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
                annotations.append(metadata)
                print(" OK")
            except Exception as e:
                print(" FAILED")
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
