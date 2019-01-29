# ARES-GT
CRISPR sgRNAs guide analysis software

## Overview

Although several online tools exist for identification of CRISPR targets (Cas9 or Cas12a/Cpf1), but usually only a limited list of genomes are available to be selected. Moreover, a limited list of candidates, based on a matrix score, is provided. I have developed ARES-GT to search CRISPR targets in DNA sequences against a reference genome provided by the user, meaning unpublished genome, unassembled, etc., can be used.

ARES-GT is licensed under the GNU Affero General Public License v3.0

## Installation

ARES-GT is a python (v2.7) script that uses datetime, sys, re and regex libraries.

## Input

  * DNA sequences in FASTA format.
  * Reference GENOME files (one file per chromosome/contig).
  * LIST with the names of reference genome files.
  
  Examples on input files are in "example" folder.
  
  ```
  If all reference genome is in FASTA format file, SplitChrom.py script will split it in files
  and will generate the list with the names of the files, ready to be used with ARES-GT:
  
  Example:
             SplitChrom.py Genome.fas
    
  ```

# Using ARES-GT

  Parameters
  
  * "f1" Precedes Name of the file with the DNA sequences in fasta format to be analysed searching for CRISPR targets.
  * "f2" Precedes Name of the file that contains the list of files with the reference genome
  * "L0" Precedes Threshold of global mismatches when number of mismatches in seed sequence in zero.
  * "L1" Precedes Threshold of global mismatches when number of mismatches in seed sequence in one.
  * "ENZ" Precedes identification of endonuclease to be used: "Cas9" or "Cas12a"
  * "NAG" Precedes "Yes" or "No", indicating is sequences with PAM "NAG" must be taken into account as possible off-targets 

```
Command line example:

ARES_GT_V1.1.py f1 Sequences.txt f2 ChromList.txt L0 5 L1 4 ENZ Cas9 NAG No 

ARES_GT_V1.1.py f1 Sequences.txt f2 ChromList.txt L0 5 L1 4 ENZ Cas12a

ARES_GT_V1.1.py f1 Genes.txt f2 Contigs.txt L0 4 L1 2 ENZ Cas9 NAG Yes
```

## Output

ARES-GT generates a tabulated text file with results.

 - 1.- A list with all the CRISPR candidate targets found in DNA sequences
``` 
ALL possible Targets found in sequences:

 SeqName	PosStart	PosEnd	Sense	Sequence(with PAM)
PDS3_EXON1_001	6	28	+	TTGGAGAAAATGGTTGTGTT TGG
PDS3_EXON1_002	7	29	+	TGGAGAAAATGGTTGTGTTT GGG
PDS3_EXON1_003	20	42	+	TGTGTTTGGGAATGTTTCTG CGG
PDS3_EXON1_004	42	64	+	GCGAATTTGCCTTATCAAAA CGG
...
```
 - 2.- A list with the optimal CRISPR candidate targets (those with no expected offtargets)
 ```
 Selected Target sequences with Unique Hits:

SeqName	PosStart	PosEnd	Sense	Sequence(with PAM)
PDS3_EXON1_006	50	72	+	GCCTTATCAAAACGGGTTTT TGG
PDS3_EXON1_011	84	106	+	TCTGGAGGTTGTGAACTAAT GGG
PDS3_EXON2_029	37	59	-	CCA AGGCCAGAGCTAGAGAACAC
PDS3_EXON2_031	95	117	-	CCT TCCGTAGTGCTCCTCGTCCT
...
```
  - 3.- Detailed information of all CRISPR candidate targets with alignment with genome location of possible offtargets
```
 <<< TARGET PDS3_EXON1_001 >>>
	PosStart	PosEnd	Sense	Sequence(with PAM)	Reverse
	6	28	+	TTGGAGAAAATGGTTGTGTT TGG	(CCA AACACAACCATTTTCTCCAA)

<Hits with NO mismatch in Seed sequence (11 nt) and LESS than 6 global mismatches.>
--When the PAM of the hit is NAG selection is LESS than 5 global mismatches:

>>>><<<<
>>HITS<<
>>>><<<<
Chrom	PosStart	PosEnd	Sense	Sequence(with PAM)	SeqIdentity	Mismatches	PAM
TAIR10_chr4	8194756	8194778	-	CCA AACACAACCATTTTCTCCAA	*** ********************	0	   NGG
TAIR10_chr1	3641599	3641621	-	CTT AACACAACCATTATATGCAG	*TT ************A*A*G**G	4	   NAG
TAIR10_chr2	4434118	4434140	-	CCA AACACAACCATTGCCTGAAA	*** ************GC**GA**	4	   NGG
TAIR10_chr2	4601042	4601064	+	TTCCAGGCAATGGTTGTGTT TGG	**CC**GC************ ***	4	   NGG
TAIR10_chr5	11006110	11006132	+	TTCCAGGCAATGGTTGTGTT TGG	**CC**GC************ ***	4	   NGG
TAIR10_chr3	6843347	6843369	-	CCC AACACAACCATTGGCCCTTA	**C ************GG*C*TT*	5	   NGG
TAIR10_chr5	4235136	4235158	-	CCT AACACAACCATCATCATCAT	**T ***********CA**AT**T	5	   NGG


<Hits with 1 mismatch in Seed sequence (11 nt) and LESS than 5 global mismatches.>
--When the PAM of the hit is NAG selection is LESS than 4 global mismatches:

>>>><<<<
>>HITS<<
>>>><<<<
Chrom	PosStart	PosEnd	Sense	Sequence(with PAM)	SeqIdentity	Mismatches	PAM
TAIR10_chr4	11335304	11335326	-	CTT AACACAATCATTTCCTCCAA	*TT *******T*****C******	2	   NAG
TAIR10_chr1	24347974	24347996	+	TTTAAGAATATGGTTGTTTT TGG	**TA****T********T** ***	4	   NGG
TAIR10_chr5	14102927	14102949	-	CCC AACACCACCATCTTCTCATA	**C *****C*****C*****AT*	4	   NGG
TAIR10_chr5	21945753	21945775	+	TTAGAGATTATGGTTTTGTT AGG	**A****TT******T**** A**	4	   NGG
```

Análisis of Results

  - List of optimal candidates (with defined parameters) is the best option when searching for CRISPR targets with low probability of offtarget mutations.
  - If no optimal candidates are found or not in the location is desired, detailed information and alignment with possible offtargets can be manually evaluated.
  - When working with gene families is expected some candidates are not selected as optimal because they match other members of the family. Those candidates can be very interesting for targeting more than one gene, which can be identify using detailed information. 


## ARES-GT algorithm

  It has been described that seed sequence (proximal sequence to PAM) is critical in the DNA binding by Cas endonucleases. In the case of Cas9 it has beed reported that seed size is around 11 nucleotides while 6-8 nucleotides in the case of Cas12a/Cpf1. Cas12a seems to be more restrictive to mismatches in the target sequence that Cas9. Additionally, distance to PAM in Cas9 and number of mismatches (and if they are grouped) have different influence in Cas9 activity. It has been also described that PAM "NAG" can be also used by Cas9 although with less efficiency.

Based on Cas9 and Cas12a sequence affinity information ARES-GT algorithm uses two steps for off-targets evaluation:

  1) All possible CRISPR target candidates are identidied in DNA sequences.
  2) Each candidate is compared with Reference Genome to evaluate possible off-targets based on:
       - Any sequence with two or more mismatches in seed sequence is discarded.
       - Any sequence with 1 mismatch in seed sequences is discarded if total mismatches are equal or more than L1 parameter.
       - Any sequence with no mismatch in seed sequences is discarded if total mismatches are equal or more than L0 parameter.



 ### Bibliography
  
-  Swarts, D. C. and M. Jinek (2018). "Cas9 versus Cas12a/Cpf1: Structure–function comparisons and implications for genome editing." Wiley Interdisciplinary Reviews: RNA 9(5): e1481.
-  Swarts, D. C., J. van der Oost and M. Jinek (2017). "Structural Basis for Guide RNA Processing and Seed-Dependent DNA Targeting by CRISPR-Cas12a." Mol Cell 66(2): 221-233 e224.
-  Sternberg, S. H., S. Redding, M. Jinek, E. C. Greene and J. A. Doudna (2014). "DNA interrogation by the CRISPR RNA-guided endonuclease Cas9." Nature.
-  Hsu, P. D., D. A. Scott, J. A. Weinstein, F. A. Ran, S. Konermann, V. Agarwala, Y. Li, E. J. Fine, X. Wu, O. Shalem, T. J. Cradick, L. A. Marraffini, G. Bao and F. Zhang (2013). "DNA targeting specificity of RNA-guided Cas9 nucleases." Nat Biotechnol 31(9): 827-832.
