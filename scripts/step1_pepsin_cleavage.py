"""
Step 1 - Fetch IsPETase sequence and predict pepsin cleavage sites
Project: Microplastic Degradation Research Project
Author:  Mamatha Gadipudi
Date:    April 2026
"""
import os, re, urllib.request
import pandas as pd
from Bio import SeqIO
from pyteomics import parser

SCRIPT_DIR   = os.path.dirname(os.path.abspath(__file__))
PROJECT_ROOT = os.path.dirname(SCRIPT_DIR)
DATA_DIR     = os.path.join(PROJECT_ROOT, "data")
RESULTS_DIR  = os.path.join(PROJECT_ROOT, "results", "step1_pepsin")
FASTA_OUT    = os.path.join(DATA_DIR,    "A0A0K8P6T7_IsPETase_wildtype.fasta")
TABLE_OUT    = os.path.join(RESULTS_DIR, "pepsin_cleavage_sites.tsv")

print("=" * 60)
print("STEP 1: Fetching IsPETase sequence from UniProt")
print("=" * 60)
ACCESSION = "A0A0K8P6T7"
URL = f"https://rest.uniprot.org/uniprotkb/{ACCESSION}.fasta"
if os.path.exists(FASTA_OUT):
    print(f"  Already downloaded: {FASTA_OUT}")
else:
    print(f"  Downloading from: {URL}")
    req = urllib.request.Request(URL, headers={"User-Agent": "Python/urllib"})
    with urllib.request.urlopen(req, timeout=30) as r:
        txt = r.read().decode("utf-8")
    with open(FASTA_OUT, "w") as f:
        f.write(txt)
    print(f"  Saved to: {FASTA_OUT}")

record = SeqIO.read(FASTA_OUT, "fasta")
wt_seq = str(record.seq)
print(f"  Protein ID     : {record.id}")
print(f"  Sequence length: {len(wt_seq)} amino acids")
print(f"  First 60 aa    : {wt_seq[:60]}")

print("\n" + "=" * 60)
print("STEP 2: Scanning for pepsin cleavage sites")
print("=" * 60)
RULE_NAME  = "pepsin ph1.3"
rule_regex = parser.expasy_rules[RULE_NAME]
print(f"  Rule   : {RULE_NAME}")
print(f"  Pattern: {rule_regex}")

sites = []
for match in re.finditer(rule_regex, wt_seq):
    p1_index    = match.start() + len(match.group()) - 1
    p1_position = p1_index + 1
    p1prime_pos = p1_position + 1
    if p1_position < 1 or p1prime_pos > len(wt_seq):
        continue
    ctx_start = max(0, p1_index - 3)
    ctx_end   = min(len(wt_seq), p1_index + 5)
    context   = wt_seq[ctx_start:ctx_end]
    sites.append({
        "cut_after_position" : p1_position,
        "P1_amino_acid"      : wt_seq[p1_index],
        "P1prime_position"   : p1prime_pos,
        "P1prime_amino_acid" : wt_seq[p1_index+1] if p1_index+1 < len(wt_seq) else "",
        "local_context"      : context,
    })

df = (pd.DataFrame(sites)
        .drop_duplicates(subset=["cut_after_position"])
        .sort_values("cut_after_position")
        .reset_index(drop=True))

df.to_csv(TABLE_OUT, sep="\t", index=False)

print(f"  Total cleavage sites found  : {len(df)}")
print(f"  Sites with Leucine (L)      : {(df['P1_amino_acid']=='L').sum()}")
print(f"  Sites with Phenylalanine (F): {(df['P1_amino_acid']=='F').sum()}")
print(f"  Saved to: {TABLE_OUT}")
print("\n  Full cleavage site table:")
print(df.to_string(index=False))
print("\n" + "=" * 60)
print("Phase 1 complete.")
print("=" * 60)
