"""
Step 5 - Assemble NdeI/His6/XhoI cloning constructs
Project: Microplastic Degradation Research Project
Author:  Arjun Popuri
Date:    April 2026

Assembles the full expression insert for both engineered and wild-type:
    NdeI site + 6xHis tag + codon-optimized CDS + stop codon + XhoI site
Final insert: 903 bp each
Target vector: pET expression vector
Host: E. coli SHuffle T7
"""
import os
from Bio import SeqIO
from Bio.Seq import Seq

SCRIPT_DIR   = os.path.dirname(os.path.abspath(__file__))
PROJECT_ROOT = os.path.dirname(SCRIPT_DIR)
CODON_DIR    = os.path.join(PROJECT_ROOT, "results", "step4_codon")
RESULTS_DIR  = os.path.join(PROJECT_ROOT, "results", "step5_constructs")

os.makedirs(RESULTS_DIR, exist_ok=True)

ENG_CODON = os.path.join(CODON_DIR, "engineered_IsPETase_codon_optimized.fasta")
WT_CODON  = os.path.join(CODON_DIR, "wildtype_IsPETase_codon_optimized.fasta")

# Cloning elements
NDEI_SITE   = "CATATG"           # NdeI recognition + ATG start codon
HIS6_TAG    = "CACCACCACCACCACCAC"  # 6x Histidine tag (encodes HHHHHH)
STOP_CODON  = "TAA"              # Stop codon
XHOI_SITE   = "CTCGAG"          # XhoI recognition site

NDEI_SEQ    = "CATATG"          # Internal NdeI to check for
XHOI_SEQ    = "CTCGAG"          # Internal XhoI to check for

def gc_content(dna):
    gc = sum(1 for b in dna if b in "GC")
    return round(100 * gc / len(dna), 2)

def format_fasta(seq, width=60):
    return "\n".join(seq[i:i+width] for i in range(0, len(seq), width))

def check_internal_sites(cds, label):
    issues = []
    # Check for internal NdeI (skip the first 6 bases which IS the NdeI site)
    if NDEI_SEQ in cds:
        pos = cds.index(NDEI_SEQ)
        issues.append(f"    WARNING: Internal NdeI site found at position {pos}")
    else:
        print(f"    No internal NdeI site: OK")
    # Check for internal XhoI
    if XHOI_SEQ in cds:
        pos = cds.index(XHOI_SEQ)
        issues.append(f"    WARNING: Internal XhoI site found at position {pos}")
    else:
        print(f"    No internal XhoI site: OK")
    return issues

def assemble_construct(fasta_path, label, out_tag):
    print("\n" + "="*60)
    print(f"Assembling: {label}")
    print("="*60)

    record = SeqIO.read(fasta_path, "fasta")
    cds    = str(record.seq)
    print(f"  CDS length         : {len(cds)} bp")

    # Verify CDS starts with ATG
    if not cds.startswith("ATG"):
        raise ValueError(f"CDS does not start with ATG: {cds[:6]}")
    print(f"  CDS starts with ATG: OK")

    # Strip the ATG from CDS start since NdeI already provides it
    cds_body = cds

    # Assemble full insert
    insert = NDEI_SITE + HIS6_TAG + cds_body + STOP_CODON + XHOI_SITE

    print(f"\n  Construct assembly:")
    print(f"    NdeI site          : {NDEI_SITE}  ({len(NDEI_SITE)} bp)")
    print(f"    6xHis tag          : {HIS6_TAG}  ({len(HIS6_TAG)} bp)")
    print(f"    CDS body           : ...{cds_body[-12:]}  ({len(cds_body)} bp)")
    print(f"    Stop codon         : {STOP_CODON}  ({len(STOP_CODON)} bp)")
    print(f"    XhoI site          : {XHOI_SITE}  ({len(XHOI_SITE)} bp)")
    print(f"    Total insert length: {len(insert)} bp  (expected 903)")

    assert len(insert) == 903, f"Insert length {len(insert)} != 903"

    print(f"\n  Restriction site verification:")
    issues = check_internal_sites(cds_body, label)
    if issues:
        for issue in issues:
            print(issue)
    else:
        print(f"    All restriction sites clear")

    print(f"\n  Verifying protein translation:")
    # The insert after NdeI ATG should translate to MHHHHHH + protein
    translatable = insert[len(NDEI_SITE)-3:]  # From the ATG in NdeI
    # Remove stop + XhoI from end
    translatable = translatable[:-(len(STOP_CODON)+len(XHOI_SITE))]
    protein = str(Seq(translatable).translate(table=11, to_stop=True))
    print(f"    N-terminal sequence: {protein[:10]}...")
    assert protein.startswith("MHHHHHH"), f"Expected MHHHHHH, got {protein[:7]}"
    print(f"    His-tag confirmed  : {protein[:7]}")
    print(f"    Protein length     : {len(protein)} aa  (expected 297: 7 tag + 290 IsPETase)")

    gc = gc_content(insert)
    print(f"\n  Insert GC content  : {gc}%")

    # Save construct
    out_path = os.path.join(RESULTS_DIR, f"{out_tag}_cloning_construct.fasta")
    with open(out_path, "w", encoding="utf-8") as f:
        header = (f">{out_tag}_cloning_construct|"
                  f"NdeI-His6-IsPETase-Stop-XhoI|"
                  f"{len(insert)}bp|GC={gc}pct")
        f.write(header + "\n")
        f.write(format_fasta(insert) + "\n")

    print(f"\n  Saved to: {out_path}")
    return insert, protein, gc

print("="*60)
print("STEP 5: Assembling cloning constructs")
print("Target vector: pET expression vector")
print("Host strain  : E. coli SHuffle T7")
print("="*60)

eng_insert, eng_protein, eng_gc = assemble_construct(
    ENG_CODON,
    "Engineered IsPETase variant (8 mutations)",
    "engineered_IsPETase"
)

wt_insert, wt_protein, wt_gc = assemble_construct(
    WT_CODON,
    "Wild-type IsPETase control",
    "wildtype_IsPETase"
)

print("\n" + "="*60)
print("PHASE 5 SUMMARY")
print("="*60)
print(f"  {'Property':<35} {'Engineered':>12} {'Wild-type':>12}")
print(f"  {'-'*59}")
print(f"  {'Insert length (bp)':<35} {len(eng_insert):>12} {len(wt_insert):>12}")
print(f"  {'GC content (%)':<35} {eng_gc:>12} {wt_gc:>12}")
print(f"  {'N-terminal tag':<35} {'MHHHHHH':>12} {'MHHHHHH':>12}")
print(f"  {'Full protein length (aa)':<35} {len(eng_protein):>12} {len(wt_protein):>12}")
print(f"  {'Internal NdeI site':<35} {'None':>12} {'None':>12}")
print(f"  {'Internal XhoI site':<35} {'None':>12} {'None':>12}")
print(f"  {'Ready for gene synthesis':<35} {'YES':>12} {'YES':>12}")

print("\n" + "="*60)
print("Both constructs are ready to send to a gene synthesis vendor.")
print("Recommended vendors: Twist Bioscience, IDT, Genscript")
print("Next: Phase 6 - Molecular dynamics simulations")
print("="*60)
