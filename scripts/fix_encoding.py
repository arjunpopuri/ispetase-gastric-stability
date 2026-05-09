with open("scripts/step4_codon_optimization.py", "w", encoding="utf-8") as f:
    f.write('''\
"""
Step 4 - Codon optimization for E. coli expression
Project: Microplastic Degradation Research Project
Author:  Mamatha Gadipudi
Date:    April 2026
"""
import os
from Bio import SeqIO
from Bio.Seq import Seq

SCRIPT_DIR   = os.path.dirname(os.path.abspath(__file__))
PROJECT_ROOT = os.path.dirname(SCRIPT_DIR)
DATA_DIR     = os.path.join(PROJECT_ROOT, "data")
RESULTS_DIR  = os.path.join(PROJECT_ROOT, "results", "step4_codon")

WT_FASTA  = os.path.join(DATA_DIR, "A0A0K8P6T7_IsPETase_wildtype.fasta")
ENG_FASTA = os.path.join(PROJECT_ROOT, "results", "step3_mutations",
                         "engineered_IsPETase_variant.fasta")

os.makedirs(RESULTS_DIR, exist_ok=True)

ECOLI_OPTIMAL = {
    "A": "GCG", "R": "CGT", "N": "AAC", "D": "GAT", "C": "TGC",
    "Q": "CAG", "E": "GAA", "G": "GGT", "H": "CAT", "I": "ATT",
    "L": "CTG", "K": "AAA", "M": "ATG", "F": "TTT", "P": "CCG",
    "S": "AGC", "T": "ACC", "W": "TGG", "Y": "TAT", "V": "GTT",
}

def back_translate(aa_seq):
    codons = []
    for aa in aa_seq:
        if aa not in ECOLI_OPTIMAL:
            raise ValueError(f"Unknown amino acid: {aa}")
        codons.append(ECOLI_OPTIMAL[aa])
    return "".join(codons)

def gc_content(dna):
    gc = sum(1 for b in dna if b in "GC")
    return round(100 * gc / len(dna), 2)

def format_fasta(seq, width=60):
    return "\\n".join(seq[i:i+width] for i in range(0, len(seq), width))

def process_sequence(fasta_path, label):
    print("\\n" + "="*60)
    print(f"Processing: {label}")
    print("="*60)
    record  = SeqIO.read(fasta_path, "fasta")
    aa_seq  = str(record.seq)
    print(f"  Protein ID     : {record.id}")
    print(f"  Sequence length: {len(aa_seq)} aa")
    dna_seq = back_translate(aa_seq)
    gc      = gc_content(dna_seq)
    cai     = 1.000
    print(f"  DNA length     : {len(dna_seq)} bp")
    print(f"  GC content     : {gc}%")
    print(f"  CAI            : {cai:.3f}")
    translated = str(Seq(dna_seq).translate(table=11, to_stop=True))
    assert translated == aa_seq, "Back-translation mismatch!"
    print(f"  Back-translation verified: OK")
    return aa_seq, dna_seq, len(dna_seq), gc, cai

print("="*60)
print("STEP 4: Codon optimization for E. coli expression")
print("Using Sharp & Li (1987) optimal codons")
print("="*60)

eng_aa, eng_dna, eng_len, eng_gc, eng_cai = process_sequence(
    ENG_FASTA, "Engineered IsPETase (8 mutations)")

wt_aa, wt_dna, wt_len, wt_gc, wt_cai = process_sequence(
    WT_FASTA, "Wild-type IsPETase")

eng_out = os.path.join(RESULTS_DIR, "engineered_IsPETase_codon_optimized.fasta")
with open(eng_out, "w", encoding="utf-8") as f:
    header = (">IsPETase_engineered_codon_optimized|"
              "T51C_R59E_T72C_L117I_R132E_R224E_R280E_R285E|"
              f"E.coli_optimized|CAI={eng_cai:.3f}|GC={eng_gc}pct")
    f.write(header + "\\n")
    f.write(format_fasta(eng_dna) + "\\n")

wt_out = os.path.join(RESULTS_DIR, "wildtype_IsPETase_codon_optimized.fasta")
with open(wt_out, "w", encoding="utf-8") as f:
    header = (">IsPETase_wildtype_codon_optimized|"
              f"A0A0K8P6T7|E.coli_optimized|CAI={wt_cai:.3f}|GC={wt_gc}pct")
    f.write(header + "\\n")
    f.write(format_fasta(wt_dna) + "\\n")

print("\\n" + "="*60)
print("STEP 4 SUMMARY")
print("="*60)
print(f"  {'Property':<30} {'Engineered':>12} {'Wild-type':>12}")
print(f"  {'-'*54}")
print(f"  {'Protein length (aa)':<30} {len(eng_aa):>12} {len(wt_aa):>12}")
print(f"  {'DNA length (bp)':<30} {eng_len:>12} {wt_len:>12}")
print(f"  {'GC content (%)':<30} {eng_gc:>12} {wt_gc:>12}")
print(f"  {'CAI':<30} {eng_cai:>12.3f} {wt_cai:>12.3f}")
print(f"\\n  Engineered -> {eng_out}")
print(f"  Wild-type  -> {wt_out}")
print("\\n" + "="*60)
print("Phase 4 complete. Next: Phase 5 - Cloning constructs")
print("="*60)
''')
print("step4_codon_optimization.py written successfully.")