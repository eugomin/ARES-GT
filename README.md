
# NEW in version 2.1

## The necessity in previous versions of having each contig of reference genome in idependent files (SplitChromV1.3.py is used in previous versions for that purpose) have been removed. Now, a unique file in FASTA format with all contigs of reference genome is used as input.


# ARES-GT (version 2.0/2.0.1/2.1)

## Overview

ARES-GT is a python script that identify CRISPR targets in query sequences (in FASTA format) and can also identify CRISPR targets matching multiple genes. Several online CRISPR tools exists with very good performance, however only a set of availables genomes can be selected (in some cases only a few, though in some tools more than 1500 can be selected) but none of them can identify sequences to be able to target multiple genes. This work usually needs manual (an slow) evaluation of alignments and offtargets.

ARES-GT run locally thus user can provide any genome as reference: unpublished data, unassembled genomes in contigs, individuals genomes for personalized treatments, etc.

ARES-GT also evaluates whether the possible CRISPR targets in any query sequence also match the other query sequences, so if a gene family is used, its algorithm will list those candidates that can target more than one members of the family.

ARES-GT is licensed under the GNU Affero General Public License v3.0

## Installation

ARES-GT is a python script that uses datetime, sys, re and regex libraries.

* Version 2.0   is a python v2.7 script.
* Version 2.0.1 has been modified for running in python 3.x versions.
* Version 2.1   use a unique fasta file as input for reference genome (run in python 3.x versions) 

## Input

  * DNA sequences in FASTA format.
  * Reference GENOME files (one file per chromosome/contig) [versions 2.0 and 2.0.1].
  * Reference GENOME file (FASTA format) [versions 2.1]
  * LIST with the names of reference genome files [versions 2.0 and 2.0.1].
  
  Examples of input files are in "example" folder.
  
  ```
  If all reference genome is in FASTA format file, SplitChrom.py script will split it in files
  and will generate the list with the names of the files, ready to be used with ARES-GT:
  
  Example:
             SplitChromV1.3.py Genome.fas
    
  ```

# Using ARES-GT (v2.1)

Parameters
  
  * "-f1" Precedes Name of the file with the DNA sequences in fasta format to be analysed searching for CRISPR targets.
  * "-f2" Precedes Name of the file that contains the reference genome in FASTA format.
  * "-L0" Precedes Threshold of global mismatches when number of mismatches in seed sequence in zero.
  * "-L1" Precedes Threshold of global mismatches when number of mismatches in seed sequence in one.
  * "-ENZ" Precedes identification of endonuclease to be used: "Cas9" or "Cas12a"
  * "-NAG" Precedes "Yes" or "No", indicating is sequences with PAM "NAG" must be taken into account as possible off-targets
  * "-OR" Preceeds "Yes" or "No". If "Yes" all CRISPR targets that do not match more than one query sequence will be discarded.

```
Command line example:

ARES_GT_V2.1.py -f1 Sequences.txt -f2 HsGRCh38.p13.fna -L0 4 -L1 3 -ENZ Cas9 -NAG No -OR No

ARES_GT_V2.1.py -f1 Sequences.fas -f2 Osativa_323_v7.0.fa -L0 4 -L1 3 -ENZ Cas12a -OR Yes

ARES_GT_V2.1.py -f1 Genes.txt -f2 OwnGenome.consensus -L0 4 -L1 2 -ENZ Cas9 -NAG Yes -OR Yes
```

# Using ARES-GT (v2.0/2.0.1)


  Parameters
  
  * "-f1" Precedes Name of the file with the DNA sequences in fasta format to be analysed searching for CRISPR targets.
  * "-f2" Precedes Name of the file that contains the list of files with the reference genome
  * "-L0" Precedes Threshold of global mismatches when number of mismatches in seed sequence in zero.
  * "-L1" Precedes Threshold of global mismatches when number of mismatches in seed sequence in one.
  * "-ENZ" Precedes identification of endonuclease to be used: "Cas9" or "Cas12a"
  * "-NAG" Precedes "Yes" or "No", indicating is sequences with PAM "NAG" must be taken into account as possible off-targets
  * "-OR" Preceeds "Yes" or "No". If "Yes" all CRISPR targets that do not match more than one query sequence will be discarded.

```
Command line example:

ARES_GT_V2.0.py -f1 Sequences.txt -f2 ChromList.txt -L0 4 -L1 3 -ENZ Cas9 -NAG No -OR No

ARES_GT_V2.0.1.py -f1 Sequences.fas -f2 ChromList.txt -L0 4 -L1 3 -ENZ Cas12a -OR Yes

ARES_GT_V2.0.py -f1 Genes.txt -f2 Contigs.txt -L0 4 -L1 2 -ENZ Cas9 -NAG Yes -OR Yes
```

# Practical examples

In the "Example V2.0" folder all files needed to test the software are available. File "Genome.txt" contains seveal contigs in fasta format, as an example of a reference genome. The first step is to extract each contig (in complete genomes it means each chromosome and any unlocated contig). The tool SplitChromV1.3.py will do the work as follow:

```
SplitChromV1.3.py Genome.txt
```

It generates one txt file for each contig and a file named "ChromosomeList.txt" with the list of all the files that have been generated. This step is recommended as contig names are reduced to 20 characters because in many genomes the name line in tha fasta format can contain many information and be too long to allow an easy display of ARES-GT results. In the case of duplicated names a number index is added with an order corresponding with the original fasta file order. In the next version that function will be automatically added to ARES-GT. (ADDED function!! See version 2.1)

Optional uses of SplitChromV1.3.py with version 2.1

	This tools generates an independent file with each contig in the reference genome file, with can be used to resolve some questions when revising genomes.
	
	1) Version 2.1 use the complete reference genome in a fasta format, however contigs/chromosome names can contains many information so name it is trimmed to the first 20 characters. If some contigs share the same first 20 characters, an index is added to the name. Using SplitChromV1.3.py the independent files names will match exactly the same names used in the ARES-GT resuts file.

	2) Prices for NGS services are getting down each month which makes easier to obtain de novo genomes of your favourite species. Anyway, a big number of contigs are provided specially when no public good reference genome is available, which is usually delivered in an unique file in fasta format. This tool easily splits that delivered file into independent files, one per contig. Then this contigs can be evaluated independently, for instance by contig-to-contig alignment and plotting. In the case of hybrid individuals, this methodology can help to identify groups of contigs that correspond with each haploid genome.



```
Practical cases
```

Case 1)
	A set of genes (in file "Genes.txt") must be evaluated and we want both specific targets for each gene and any target matching several genes. The selected endonuclease is Cas9 and offtargets with "NAG" PAM must be taken into account. User decide that offtarget threshold will be 5 total mismatches when no mismatches are found in seed sequences and 3 total mismatches when 1 mismatch is found in seed sequence:

```
ARES_GT_V2.0.py -f1 Genes.txt -f2 ChromosomeList.txt -L0 5 -L1 3 -ENZ Cas9 -NAG Yes -OR No
```

Case 2)
	The same parameters than "case 1" but we only want to evaluate possible targets that match several of the sequences in Genes.txt file, discarding targets only matching one sequence:
	
```
ARES_GT_V2.0.py -f1 Genes.txt -f2 ChromosomeList.txt -L0 5 -L1 3 -ENZ Cas9 -NAG Yes -OR Yes
```

Case 3)
	The selected endonuclease is Cas12a and total mismatches threshold in the case of no mismatches in seed sequence will be 4 and 3 in the case of one mismatch in seed sequence. We only want possible targets matching several query sequences.
	
```
ARES_GT_V2.0.py -f1 Genes.txt -f2 ChromosomeList.txt -L0 4 -L1 3 -ENZ Cas12a -OR Yes
```

## Output (new in v2.0/2.0.1)

ARES-GT generates a tabulated text file with results. The output format is similar to previous versions (see below in v1.2 output explanation). Some repetitive information have been also eliminated to make easier the results file interpretation.

```
[[[ TARGET KYO_CBF1_017 ]]]
	PosStart	PosEnd	Sense	Sequence(with PAM)	Reverse
	510	532	+	TGACGAACTCCTCTGTAAAT TGG	(CCA ATTTACAGAGGAGTTCGTCA)
>>>><<<<
>>HITS<<  NO mismatch in Seed sequence (11 nt)
>>>><<<<

Chrom	PosStart	PosEnd	Sense	Sequence(with PAM)	SeqIdentity	Mismatches	PAM
TAIR10_chr4	13015920	13015942	+	TGACGAACTCCTCTGTAAAT TGG	******************** ***	0	   NGG
TAIR10_chr4	13022405	13022427	+	TGACGAACTCCTCTGTAAAT TGG	******************** ***	0	   NGG
TAIR10_chr5	21117612	21117634	+	TGACGAACTCCTCTGTAAAT CGG	******************** C**	0	   NGG

>>>><<<<
>>HITS<<  1 mismatch in Seed sequence (11 nt)
>>>><<<<

Chrom	PosStart	PosEnd	Sense	Sequence(with PAM)	SeqIdentity	Mismatches	PAM
TAIR10_chr4	13018837	13018859	+	CGACGAACTCCTCTGTATAT TGG	C****************T** ***	2	   NGG
```

New information have been added for listing those CRISPR sequences that can targets more than one query sequence. At the end of each line the number of possible off-targets (taking into account the multiple targets) is indicated. 

```
Example of candidate matching several sequences:

KYO_CBF1_014 = KYO_CBF2_081	AGCACGAGCTGCCATCTCAG  NGG	Number of Possible OffTargets = 1
KYO_CBF1_015 = KYO_CBF2_082	GAGCTGCCATCTCAGCGGTT  NGG	Number of Possible OffTargets = 0
KYO_CBF1_017 = KYO_CBF2_084 = KYO_CBF4_211	TGACGAACTCCTCTGTAAAT  NGG	Number of Possible OffTargets = 1
KYO_CBF1_018 = KYO_CBF2_085	GACGAACTCCTCTGTAAATT  NGG	Number of Possible OffTargets = 2
KYO_CBF1_019 = KYO_CBF3_153	CGAAACTTCTTACGACCCGC  NGG	Number of Possible OffTargets = 0
... ...
```

## In previous Versions:

# ARES-GT (version 1.2.1)
CRISPR sgRNAs guide analysis software

## Overview

Although several online tools exist for identification of CRISPR targets (Cas9 or Cas12a/Cpf1), usually only a limited list of genomes are available to be selected. Moreover, a limited list of candidates, based on a matrix score, is provided. I have developed ARES-GT to search CRISPR targets in DNA sequences against a reference genome provided by the user, meaning unpublished genome, unassembled, etc., can be used.

ARES-GT is licensed under the GNU Affero General Public License v3.0

## Installation

ARES-GT is a python (v2.7) script that uses datetime, sys, re and regex libraries.

## Input

  * DNA sequences in FASTA format.
  * Reference GENOME files (one file per chromosome/contig).
  * LIST with the names of reference genome files.
  
  Examples of input files are in "example" folder.
  
  ```
  If all reference genome is in FASTA format file, SplitChrom.py script will split it in files
  and will generate the list with the names of the files, ready to be used with ARES-GT:
  
  Example:
             SplitChrom.py Genome.fas
    
  ```

# Using ARES-GT (v1.2.1)

  Parameters
  
  * "-f1" Precedes Name of the file with the DNA sequences in fasta format to be analysed searching for CRISPR targets.
  * "-f2" Precedes Name of the file that contains the list of files with the reference genome
  * "-L0" Precedes Threshold of global mismatches when number of mismatches in seed sequence in zero.
  * "-L1" Precedes Threshold of global mismatches when number of mismatches in seed sequence in one.
  * "-ENZ" Precedes identification of endonuclease to be used: "Cas9" or "Cas12a"
  * "-NAG" Precedes "Yes" or "No", indicating is sequences with PAM "NAG" must be taken into account as possible off-targets 

```
Command line example:

ARES_GT_V1.1.py -f1 Sequences.txt -f2 ChromList.txt -L0 5 -L1 4 -ENZ Cas9 -NAG No 

ARES_GT_V1.1.py -f1 Sequences.txt -f2 ChromList.txt -L0 5 -L1 4 -ENZ Cas12a

ARES_GT_V1.1.py -f1 Genes.txt -f2 Contigs.txt -L0 4 -L1 2 -ENZ Cas9 -NAG Yes
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

  1) All possible CRISPR target candidates are identidied in query DNA sequences.
  2) Each candidate is compared with Reference Genome to evaluate possible off-targets based on:
       - Any sequence with two or more mismatches in seed sequence is discarded.
       - Any sequence with 1 mismatch in seed sequences is discarded if the total number  mismatches is more than L1 parameter.
       - Any sequence with no mismatch in seed sequences is discarded if the total number of mismatches is more than L0 parameter.



 ### CITATION
  
If you use ARES-GT please cite:

Minguet, E. G. (2020). "Ares-GT: Design of guide RNAs targeting multiple genes for CRISPR-Cas experiments." PLOS ONE 15(10): e0241001.
https://doi.org/10.1371/journal.pone.0241001
