code = '''\
"""
Step 2 - Download PDB 5XJH and compute surface exposure (SASA)
Project: Microplastic Degradation Research Project
Author:  Arjun Popuri
Date:    April 2026
"""
import os, urllib.request, warnings
import numpy as np
import pandas as pd
from Bio.PDB import PDBParser
from Bio.PDB.SASA import ShrakeRupley
from Bio.PDB.Polypeptide import protein_letters_3to1

warnings.filterwarnings("ignore")

SCRIPT_DIR   = os.path.dirname(os.path.abspath(__file__))
PROJECT_ROOT = os.path.dirname(SCRIPT_DIR)
DATA_DIR     = os.path.join(PROJECT_ROOT, "data")
RESULTS_DIR  = os.path.join(PROJECT_ROOT, "results", "step2_structure")

PDB_OUT      = os.path.join(DATA_DIR,    "5XJH_IsPETase.pdb")
ANNOT_OUT    = os.path.join(RESULTS_DIR, "residue_sasa_annotations.tsv")
CHARGE_OUT   = os.path.join(RESULTS_DIR, "charge_reversal_candidates.tsv")
DISULF_OUT   = os.path.join(RESULTS_DIR, "disulfide_candidates.tsv")

MAX_ASA = {
    "ALA":129,"ARG":274,"ASN":195,"ASP":193,"CYS":167,
    "GLN":225,"GLU":223,"GLY":104,"HIS":224,"ILE":197,
    "LEU":201,"LYS":236,"MET":224,"PHE":240,"PRO":159,
    "SER":155,"THR":172,"TRP":285,"TYR":263,"VAL":174,
}

CATALYTIC = [160, 206, 237]

print("=" * 60)
print("STEP 1: Downloading PDB structure 5XJH")
print("=" * 60)

if os.path.exists(PDB_OUT):
    print(f"  Already downloaded: {PDB_OUT}")
else:
    url = "https://files.rcsb.org/download/5XJH.pdb"
    print(f"  Downloading from: {url}")
    urllib.request.urlretrieve(url, PDB_OUT)
    print(f"  Saved to: {PDB_OUT} ({os.path.getsize(PDB_OUT):,} bytes)")

parser = PDBParser(QUIET=True)
structure = parser.get_structure("5XJH", PDB_OUT)
model   = structure[0]
chain_A = model["A"]

std_residues = [r for r in chain_A if r.id[0] == " " and r.resname in MAX_ASA]
print(f"  Standard residues in chain A: {len(std_residues)}")

print("\\n" + "=" * 60)
print("STEP 2: Computing SASA for every residue")
print("=" * 60)

sr = ShrakeRupley()
sr.compute(structure, level="R")
print("  SASA computation complete.")

cat_ca = []
for rn in CATALYTIC:
    res = chain_A[(" ", rn, " ")]
    cat_ca.append(res["CA"].get_vector().get_array())
print(f"  Catalytic triad located: Ser{CATALYTIC[0]}, Asp{CATALYTIC[1]}, His{CATALYTIC[2]}")

rows = []
for res in std_residues:
    resseq  = res.id[1]
    resname = res.resname
    one_let = protein_letters_3to1.get(resname, "X")
    asa     = res.sasa
    rel_asa = asa / MAX_ASA[resname]
    if "CA" in res:
        ca       = res["CA"].get_vector().get_array()
        min_dist = float(min(np.linalg.norm(ca - c) for c in cat_ca))
    else:
        min_dist = 999.0
    rows.append({
        "resseq"                   : resseq,
        "resname"                  : resname,
        "one_letter"               : one_let,
        "ASA"                      : round(asa, 2),
        "relative_ASA"             : round(rel_asa, 4),
        "min_dist_to_catalytic"    : round(min_dist, 2),
        "is_surface"               : rel_asa >= 0.25,
        "near_active_site"         : min_dist < 8.0,
    })

df_ann = pd.DataFrame(rows)
df_ann.to_csv(ANNOT_OUT, sep="\\t", index=False)
print(f"  Annotated {len(df_ann)} residues -> {ANNOT_OUT}")
print(f"  Surface-exposed residues (rel_ASA >= 0.25): {df_ann.is_surface.sum()}")

print("\\n" + "=" * 60)
print("STEP 3: Identifying charge-reversal candidates")
print("=" * 60)

df_charge = df_ann[
    df_ann["is_surface"] &
    ~df_ann["near_active_site"] &
    df_ann["one_letter"].isin(["R", "K"])
].copy()
df_charge["proposed_mutation"] = df_charge["one_letter"] + df_charge["resseq"].astype(str) + "E"
df_charge["delta_charge"]      = -2
df_charge.to_csv(CHARGE_OUT, sep="\\t", index=False)
print(f"  Charge-reversal candidates found: {len(df_charge)}")
print(df_charge[["resseq","one_letter","relative_ASA","min_dist_to_catalytic","proposed_mutation"]].to_string(index=False))

print("\\n" + "=" * 60)
print("STEP 4: Identifying disulfide bond candidates")
print("=" * 60)

from itertools import combinations

elig = df_ann[
    df_ann["is_surface"] &
    ~df_ann["near_active_site"] &
    ~df_ann["one_letter"].isin(["G","P","C"])
].reset_index(drop=True)

coord_map = {}
for _, row in elig.iterrows():
    rseq = int(row["resseq"])
    try:
        res = chain_A[(" ", rseq, " ")]
        if "CA" in res and "CB" in res:
            coord_map[rseq] = (
                res["CA"].get_vector().get_array(),
                res["CB"].get_vector().get_array()
            )
    except KeyError:
        pass

elig = elig[elig["resseq"].isin(coord_map)].reset_index(drop=True)

pair_rows = []
for i, j in combinations(range(len(elig)), 2):
    r1 = elig.iloc[i]; r2 = elig.iloc[j]
    rseq1 = int(r1["resseq"]); rseq2 = int(r2["resseq"])
    ca1, cb1 = coord_map[rseq1]
    ca2, cb2 = coord_map[rseq2]
    cb_dist  = float(np.linalg.norm(cb1 - cb2))
    ca_dist  = float(np.linalg.norm(ca1 - ca2))
    if not (3.5 <= cb_dist <= 5.5):
        continue
    min_triad = min(r1["min_dist_to_catalytic"], r2["min_dist_to_catalytic"])
    pair_rows.append({
        "resseq1"         : rseq1,
        "one_letter1"     : r1["one_letter"],
        "resseq2"         : rseq2,
        "one_letter2"     : r2["one_letter"],
        "CB_CB_dist_A"    : round(cb_dist, 3),
        "CA_CA_dist_A"    : round(ca_dist, 3),
        "avg_rel_ASA"     : round((r1["relative_ASA"]+r2["relative_ASA"])/2, 4),
        "min_dist_to_triad": round(min_triad, 2),
        "proposed_mutation": f"{r1['one_letter']}{rseq1}C + {r2['one_letter']}{rseq2}C",
    })

df_ds = (pd.DataFrame(pair_rows)
           .sort_values("CB_CB_dist_A")
           .reset_index(drop=True))
df_ds.to_csv(DISULF_OUT, sep="\\t", index=False)
print(f"  Disulfide candidate pairs found: {len(df_ds)}")
print(df_ds[["resseq1","resseq2","CB_CB_dist_A","avg_rel_ASA","min_dist_to_triad","proposed_mutation"]].to_string(index=False))

print("\\n" + "=" * 60)
print("Phase 2 complete.")
print(f"  Annotations -> {ANNOT_OUT}")
print(f"  Charge candidates -> {CHARGE_OUT}")
print(f"  Disulfide candidates -> {DISULF_OUT}")
print("=" * 60)
'''

with open("scripts/step2_structure_analysis.py", "w") as f:
    f.write(code)
print("step2_structure_analysis.py written successfully.")
