"""
Step 3 - Design the 8-mutation engineered IsPETase variant
Project: Microplastic Degradation Research Project
Author:  Mamatha Gadipudi
Date:    April 2026

Mutations applied:
  T51C  - disulfide bond (pairs with T72C)
  T72C  - disulfide bond (pairs with T51C)
  R59E  - charge reversal (surface Arg, high exposure)
  L117I - pepsin resistance (only surface-exposed high-risk site)
  R132E - charge reversal (surface Arg)
  R224E - charge reversal (surface Arg)
  R280E - charge reversal (surface Arg, sits in subsite IIc)
  R285E - charge reversal (surface Arg)
"""
import os
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

SCRIPT_DIR   = os.path.dirname(os.path.abspath(__file__))
PROJECT_ROOT = os.path.dirname(SCRIPT_DIR)
DATA_DIR     = os.path.join(PROJECT_ROOT, "data")
RESULTS_DIR  = os.path.join(PROJECT_ROOT, "results", "step3_mutations")

FASTA_WT  = os.path.join(DATA_DIR, "A0A0K8P6T7_IsPETase_wildtype.fasta")
FASTA_ENG = os.path.join(RESULTS_DIR, "engineered_IsPETase_variant.fasta")
TABLE_OUT = os.path.join(RESULTS_DIR, "mutation_summary.tsv")

MUTATIONS = [
    (51,  "T", "C", "disulfide",      "Forms S-S bond with T72C, Cb-Cb = 4.50 A"),
    (72,  "T", "C", "disulfide",      "Partner of T51C disulfide bond"),
    (59,  "R", "E", "charge-reversal","Highest surface exposure rel_ASA=0.68, 26.0 A from triad"),
    (117, "L", "I", "pepsin",         "Only surface-exposed high-risk pepsin cleavage site"),
    (132, "R", "E", "charge-reversal","Surface Arg, rel_ASA=0.40, 22.3 A from triad"),
    (224, "R", "E", "charge-reversal","Surface Arg, rel_ASA=0.54, 22.7 A from triad"),
    (280, "R", "E", "charge-reversal","Surface Arg in subsite IIc, rel_ASA=0.48, 15.0 A from triad"),
    (285, "R", "E", "charge-reversal","Surface Arg, rel_ASA=0.40, 15.0 A from triad"),
]

CATALYTIC_TRIAD = {160: "S", 206: "D", 237: "H"}

print("=" * 60)
print("STEP 1: Loading wild-type IsPETase sequence")
print("=" * 60)

record = SeqIO.read(FASTA_WT, "fasta")
wt_seq = str(record.seq)
print(f"  Protein ID     : {record.id}")
print(f"  Sequence length: {len(wt_seq)} aa")
assert len(wt_seq) == 290, f"Expected 290 aa, got {len(wt_seq)}"
print(f"  Length confirmed: 290 aa")

seq_map = {i+1: aa for i, aa in enumerate(wt_seq)}

print("\n" + "=" * 60)
print("STEP 2: Verifying WT residues at each mutation position")
print("=" * 60)

all_match = True
for pos, wt_aa, new_aa, category, rationale in MUTATIONS:
    actual = seq_map[pos]
    match  = actual == wt_aa
    status = "OK" if match else f"MISMATCH - expected {wt_aa}"
    if not match:
        all_match = False
    print(f"  Pos {pos:>3}: {actual} -> {new_aa}  ({category})  {status}")

if not all_match:
    raise ValueError("One or more WT residues did not match. Check mutation list.")
print("  All WT residues verified.")

print("\n" + "=" * 60)
print("STEP 3: Applying mutations")
print("=" * 60)

eng_list = list(wt_seq)
for pos, wt_aa, new_aa, category, rationale in MUTATIONS:
    eng_list[pos - 1] = new_aa
eng_seq = "".join(eng_list)

n_diffs = sum(a != b for a, b in zip(wt_seq, eng_seq))
print(f"  Mutations applied: {n_diffs} (expected 8)")

print("\n" + "=" * 60)
print("STEP 4: Verifying catalytic triad is untouched")
print("=" * 60)

for pos, expected_aa in CATALYTIC_TRIAD.items():
    actual = eng_seq[pos - 1]
    status = "OK" if actual == expected_aa else "CHANGED - PROBLEM"
    print(f"  Catalytic residue {pos} ({expected_aa}): {actual}  {status}")
    assert actual == expected_aa, f"Catalytic residue {pos} was altered!"
print("  Catalytic triad intact.")

print("\n" + "=" * 60)
print("STEP 5: Side-by-side mutation comparison")
print("=" * 60)

rows = []
print(f"  {'Pos':>4}  {'WT':>3}  {'Eng':>3}  {'Category':<16}  Rationale")
print("  " + "-" * 80)
for pos, wt_aa, new_aa, category, rationale in MUTATIONS:
    print(f"  {pos:>4}  {wt_aa:>3}  {new_aa:>3}  {category:<16}  {rationale}")
    rows.append({
        "position"        : pos,
        "WT_residue"      : wt_aa,
        "engineered_residue": new_aa,
        "category"        : category,
        "rationale"       : rationale,
    })

print("\n" + "=" * 60)
print("STEP 6: Net charge analysis at pH 1.2")
print("=" * 60)

charge_reversals = [(pos, wt, new) for pos, wt, new, cat, _ in MUTATIONS if cat == "charge-reversal"]
wt_charge_contrib  = len(charge_reversals) * 1   # each Arg = +1 at pH 1.2
eng_charge_contrib = len(charge_reversals) * (-1) # each Glu = -1 at pH 1.2
net_shift = eng_charge_contrib - wt_charge_contrib

print(f"  Charge-reversal mutations: {len(charge_reversals)}")
print(f"  Each Arg->Glu swap removes +2 charge units at pH 1.2")
print(f"  Total charge shift: {net_shift * -1} units removed")
print(f"  Estimated WT net charge at pH 1.2 : +22e")
print(f"  Estimated ENG net charge at pH 1.2: +17e")

print("\n" + "=" * 60)
print("STEP 7: Saving engineered FASTA and mutation table")
print("=" * 60)

eng_id = ("A0A0K8P6T7_IsPETase|engineered|"
          "T51C_R59E_T72C_L117I_R132E_R224E_R280E_R285E")
eng_desc = "Engineered IsPETase variant with 8 targeted mutations for gastric stability"

with open(FASTA_ENG, "w") as f:
    f.write(f">{eng_id} {eng_desc}\n")
    for i in range(0, len(eng_seq), 60):
        f.write(eng_seq[i:i+60] + "\n")

df_mut = pd.DataFrame(rows)
df_mut.to_csv(TABLE_OUT, sep="\t", index=False)

print(f"  Engineered FASTA -> {FASTA_ENG}")
print(f"  Mutation table   -> {TABLE_OUT}")
print(f"  Sequence length  : {len(eng_seq)} aa (expected 290)")

print("\n" + "=" * 60)
print("Phase 3 complete.")
print("=" * 60)
