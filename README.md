\# Computational Engineering of a Gastric-Stable IsPETase Variant

\### Oral Enzymatic Therapeutic for PET Microplastic Degradation



\*\*Researcher:\*\* Arjun Popuri  

\*\*Date:\*\* April 2026  

\*\*Location:\*\* Dallas, Texas  



\---



\## Overview



This repository contains the full computational pipeline for engineering a 

gastric acid-resistant variant of \*Ideonella sakaiensis\* PETase (IsPETase) — 

a bacterial enzyme that degrades PET plastic. The goal is to design an enzyme 

variant that survives the human stomach (pH 1.2, pepsin) and degrades ingested 

PET microplastics in the gastrointestinal tract.



Wild-type IsPETase is rapidly inactivated under gastric conditions due to two 

failure mechanisms: electrostatic destabilization from protonated surface 

Arginine residues at pH 1.2, and proteolytic degradation by pepsin at 

surface-exposed cleavage sites. This study addresses both simultaneously.



\---



\## Hypothesis



An engineered IsPETase variant incorporating five Arg→Glu charge-reversal 

mutations (R59E, R132E, R224E, R280E, R285E), an engineered T51C/T72C 

disulfide bond, and a conservative L117I pepsin resistance substitution will 

demonstrate significantly enhanced residual PET-hydrolyzing activity following 

simulated gastric fluid incubation compared to wild-type.



\---



\## The 8 Mutations



| Position | WT | Mut | Category | Rationale |

|---------|-----|-----|----------|-----------|

| 51 | T | C | Disulfide | Forms S-S bond with T72C, Cb-Cb = 4.50 A |

| 72 | T | C | Disulfide | Partner of T51C |

| 59 | R | E | Charge reversal | Highest surface exposure rel\_ASA = 0.68 |

| 117 | L | I | Pepsin resistance | Only surface-exposed high-risk pepsin site |

| 132 | R | E | Charge reversal | rel\_ASA = 0.40, 22.3 A from catalytic triad |

| 224 | R | E | Charge reversal | rel\_ASA = 0.54, 22.7 A from catalytic triad |

| 280 | R | E | Charge reversal | Subsite IIc, rel\_ASA = 0.48 |

| 285 | R | E | Charge reversal | rel\_ASA = 0.40 |



\---



\## Computational Pipeline



\### Phase 1 — Pepsin Cleavage Site Mapping

\- Retrieved wild-type IsPETase sequence from UniProt (A0A0K8P6T7, 290 aa)

\- Applied ExPASy pepsin pH 1.3 cleavage rule via pyteomics

\- Identified 35 cleavage sites; L117 confirmed as only actionable surface site

\- Output: `results/step1\_pepsin/pepsin\_cleavage\_sites.tsv`



\### Phase 2 — Structural SASA Analysis

\- Downloaded crystal structure PDB 5XJH (Joo et al. 2018, 1.5 A resolution)

\- Computed per-residue SASA using ShrakeRupley algorithm (BioPython 1.87)

\- Identified 7 charge-reversal candidates and 16 disulfide bond candidates

\- Catalytic triad confirmed: Ser160 / Asp206 / His237

\- Output: `results/step2\_structure/`



\### Phase 3 — Variant Design and Verification

\- Applied 8 mutations with pre-mutation WT residue verification

\- Post-mutation catalytic triad integrity check — all passed

\- Net surface charge reduced from +22e to +17e at pH 1.2

\- Output: `results/step3\_mutations/engineered\_IsPETase\_variant.fasta`



\### Phase 4 — Codon Optimization

\- Back-translated both sequences using Sharp and Li (1987) E. coli optimal codons

\- CAI = 1.000 for both variants (theoretical maximum)

\- GC content: 56.67% engineered, 57.47% wild-type

\- Output: `results/step4\_codon/`



\### Phase 5 — Cloning Construct Assembly

\- Assembled NdeI / 6xHis / CDS / Stop / XhoI expression constructs

\- Total insert: 903 bp each

\- Verified: no internal NdeI or XhoI sites, correct MHHHHHH N-terminal tag

\- Expression host: E. coli SHuffle T7 (required for disulfide bond formation)

\- Output: `results/step5\_constructs/`



\---



\## Key Results



| Metric | Wild-type | Engineered |

|--------|-----------|------------|

| Net surface charge at pH 1.2 | +22e | +17e |

| Charge shift | — | -10 units |

| Codon Adaptation Index | 1.000 | 1.000 |

| Insert length | 903 bp | 903 bp |

| Catalytic triad intact | Yes | Yes |

| Internal restriction sites | None | None |

| Ready for gene synthesis | Yes | Yes |



\---



\## Repository Structure

