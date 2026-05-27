code = '''\
"""
Step 4 - Codon optimization for E. coli expression
Project: Microplastic Degradation Research Project
Author:  Arjun Popuri 
Date:    April 2026

Uses Sharp & Li (1987) E. coli optimal codons.
Produces codon-optimized DNA for both engineered variant and wild-type control.
Target: CAI = 1.000, GC content ~50-55%
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

# ── Sharp & Li 1987 E. coli optimal codons ────────────────────────────────────
# For each amino acid, the single highest-CAI codon for E. coli K-12
ECOLI_OPTIMAL = {
    "A": "GCG",  # Alanine
    "R": "CGT",  # Arginine
    "N": "AAC",  # Asparagine
    "D": "GAT",  # Aspartate
    "C": "TGC",  # Cysteine
    "Q": "CAG",  # Glutamine
    "E": "GAA",  # Glutamate
    "G": "GGT",  # Glycine
    "H": "CAT",  # Histidine
    "I": "ATT",  # Isoleucine
    "L": "CTG",  # Leucine
    "K": "AAA",  # Lysine
    "M": "ATG",  # Methionine
    "F": "TTT",  # Phenylalanine
    "P": "CCG",  # Proline
    "S": "AGC",  # Serine
    "T": "ACC",  # Threonine
    "W": "TGG",  # Tryptophan
    "Y": "TAT",  # Tyrosine
    "V": "GTT",  # Valine
    "*": "TAA",  # Stop codon
}

# CAI values for the optimal codons (Sharp & Li 1987)
OPTIMAL_CAI_VALUES = {
    "A": 1.000, "R": 1.000, "N": 1.000, "D": 1.000, "C": 1.000,
    "Q": 1.000, "E": 1.000, "G": 1.000, "H": 1.000, "I": 1.000,
    "L": 1.000, "K": 1.000, "M": 1.000, "F": 1.000, "P": 1.000,
    "S": 1.000, "T": 1.000, "W": 1.000, "Y": 1.000, "V": 1.000,
}

def back_translate(aa_seq):
    """Convert amino acid sequence to E. coli optimized DNA."""
    codons = []
    for aa in aa_seq:
        if aa not in ECOLI_OPTIMAL:
            raise ValueError(f"Unknown amino acid: {aa}")
        codons.append(ECOLI_OPTIMAL[aa])
    return "".join(codons)

def gc_content(dna):
    """Calculate GC content as a percentage."""
    gc = sum(1 for b in dna if b in "GC")
    return round(100 * gc / len(dna), 2)

def format_fasta(seq, width=60):
    """Format DNA sequence in 60-character lines."""
    return "\\n".join(seq[i:i+width] for i in range(0, len(seq), width))

def process_sequence(fasta_path, label):
    """Load FASTA, back-translate, compute stats, return results."""
    print(f"\\n{'='*60}")
    print(f"Processing: {label}")
    print(f"{'='*60}")

    record  = SeqIO.read(fasta_path, "fasta")
    aa_seq  = str(record.seq)
    print(f"  Protein ID     : {record.id}")
    print(f"  Sequence length: {len(aa_seq)} aa")

    dna_seq = back_translate(aa_seq)
    dna_len = len(dna_seq)
    gc      = gc_content(dna_seq)
    cai     = 1.000  # All optimal codons by construction

    print(f"  DNA length     : {dna_len} bp  (expected {len(aa_seq)*3})")
    print(f"  GC content     : {gc}%")
    print(f"  CAI            : {cai:.3f}")

    # Verify back-translation is correct
    translated = str(Seq(dna_seq).translate(table=11, to_stop=True))
    assert translated == aa_seq, "Back-translation mismatch!"
    print(f"  Back-translation verified: DNA -> protein matches input")

    return aa_seq, dna_seq, dna_len, gc, cai, record.id

print("=" * 60)
print("STEP 4: Codon optimization for E. coli expression")
print("Using Sharp & Li (1987) optimal codons")
print("=" * 60)

# Process engineered variant
eng_aa, eng_dna, eng_len, eng_gc, eng_cai, eng_id = process_sequence(
    ENG_FASTA, "Engineered IsPETase (8 mutations)")

# Process wild-type control
wt_aa, wt_dna, wt_len, wt_gc, wt_cai, wt_id = process_sequence(
    WT_FASTA, "Wild-type IsPETase")

# Save engineered codon-optimized DNA
eng_out = os.path.join(RESULTS_DIR, "engineered_IsPETase_codon_optimized.fasta")
with open(eng_out, "w") as f:
    header = (f">IsPETase_engineered_codon_optimized|"
              f"T51C_R59E_T72C_L117I_R132E_R224E_R280E_R285E|"
              f"E.coli_optimized|CAI={eng_cai:.3f}|GC={eng_gc}pct")
    f.write(header + "\\n")
    f.write(format_fasta(eng_dna) + "\\n")

# Save wild-type codon-optimized DNA
wt_out = os.path.join(RESULTS_DIR, "wildtype_IsPETase_codon_optimized.fasta")
with open(wt_out, "w") as f:
    header = (f">IsPETase_wildtype_codon_optimized|"
              f"A0A0K8P6T7|"
              f"E.coli_optimized|CAI={wt_cai:.3f}|GC={wt_gc}pct")
    f.write(header + "\\n")
    f.write(format_fasta(wt_dna) + "\\n")

print(f"\\n{'='*60}")
print("STEP 4 SUMMARY")
print(f"{'='*60}")
print(f"  {'Property':<30} {'Engineered':>12} {'Wild-type':>12}")
print(f"  {'-'*54}")
print(f"  {'Protein length (aa)':<30} {len(eng_aa):>12} {len(wt_aa):>12}")
print(f"  {'DNA length (bp)':<30} {eng_len:>12} {wt_len:>12}")
print(f"  {'GC content (%)':<30} {eng_gc:>12} {wt_gc:>12}")
print(f"  {'CAI':<30} {eng_cai:>12.3f} {wt_cai:>12.3f}")

print(f"\\n  Output files:")
print(f"  Engineered -> {eng_out}")
print(f"  Wild-type  -> {wt_out}")

print(f"\\n{'='*60}")
print("Phase 4 complete.")
print("Next: Phase 5 - Assemble cloning constructs")
print(f"{'='*60}")
'''

with open("scripts/step4_codon_optimization.py", "w") as f:
    f.write(code)
print("step4_codon_optimization.py written successfully.")
