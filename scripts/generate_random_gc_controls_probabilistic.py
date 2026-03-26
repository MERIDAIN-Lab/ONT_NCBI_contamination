#!/usr/bin/env python3

import random

input_fasta = "ONT_and_Porechop_adapters_and_barcodes_CLEANED.fasta"
output_fasta = "random_gc_matched_controls.fasta"
replicates = 10

def read_fasta(path):
    records = []
    with open(path) as f:
        header = None
        seq = []
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    records.append((header, "".join(seq).upper()))
                header = line[1:]
                seq = []
            else:
                seq.append(line)
        if header is not None:
            records.append((header, "".join(seq).upper()))
    return records

def gc_fraction(seq):
    gc = seq.count("G") + seq.count("C")
    return gc / len(seq)

def random_gc_matched_sequence(length, gc):
    at = 1 - gc
    bases = ["A", "T", "G", "C"]
    probs = [at / 2, at / 2, gc / 2, gc / 2]
    return "".join(random.choices(bases, weights=probs, k=length))

records = read_fasta(input_fasta)

with open(output_fasta, "w") as out:
    for header, seq in records:
        seq = seq.replace("U", "T")
        seq = "".join(base for base in seq if base in "ACGT")
        if not seq:
            continue

        gc = gc_fraction(seq)

        for i in range(replicates):
            random_seq = random_gc_matched_sequence(len(seq), gc)
            out.write(f">{header}_gc_random_{i+1}\n{random_seq}\n")

print("Output written to:", output_fasta)