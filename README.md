# ARES-GT
CRISPR sgRNAs guide analysis software

ARES-GT is a command line Python software for identification of CRISPR target in DNA sequences (Cas9 or Cas12a/Cpf1).
It uses several python libraries: datetime, sys, re, regex.

Three input files are needed:

  First: DNA sequences in fasta format in a text file.
  
  Second: Reference sequences (Genome) in text files (one file per chromosome/contig). If all reference genome is in a unique file (in fasta format) the python script "SplitGenome.py" will split the file to generate one file per contig and also a file with the list of all new files (used to indicate ARES-GT which files contains reference genome).
  
    >>> Spliting chromosomes/contigs in files helps to reduce the memory use as in some species genome size can be really big.
    
   Third: A text file that contains a list of all the files that contains genome chromosomes/contigs.

ARES-GT algorithm

  First, DNA sequences are analysed and all positions with a PAM is selected as a possible CRISPR target candidate. Then evaluation of possible offtargets is performed.

  It has been described that seed sequence (proximal sequence to PAM) is critical in the DNA binding by Cas endonucleases. In the case of Cas9 it has beed reported that seed size is around 11 nucleotides while 6-8 nucleotides in the case of Cas12a/Cpf1. Cas12a seems to be more restrictive to mismatches in the target sequence that Cas9. Additionally, distance to PAM in Cas9 and number of mismatches (and if they are grouped) have different influence in Cas9 activity. It has been also described that PAM "NAG" can be also used by Cas9 although with less efficiency.
  
  Bibliography
  
  Swarts, D. C. and M. Jinek (2018). "Cas9 versus Cas12a/Cpf1: Structure–function comparisons and implications for genome editing." Wiley Interdisciplinary Reviews: RNA 9(5): e1481.
  Swarts, D. C., J. van der Oost and M. Jinek (2017). "Structural Basis for Guide RNA Processing and Seed-Dependent DNA Targeting by CRISPR-Cas12a." Mol Cell 66(2): 221-233 e224.
  Sternberg, S. H., S. Redding, M. Jinek, E. C. Greene and J. A. Doudna (2014). "DNA interrogation by the CRISPR RNA-guided endonuclease Cas9." Nature.
  Hsu, P. D., D. A. Scott, J. A. Weinstein, F. A. Ran, S. Konermann, V. Agarwala, Y. Li, E. J. Fine, X. Wu, O. Shalem, T. J. Cradick, L. A. Marraffini, G. Bao and F. Zhang (2013). "DNA targeting specificity of RNA-guided Cas9 nucleases." Nat Biotechnol 31(9): 827-832.

Based on Cas9 and Cas12a sequence affinity information ARES-GT algorithm uses two steps for off-targets evaluation:

  1) Any sequence with two or more mismatches in seed sequence is discarded.
  2) Any sequence with less mismatches (in all 20 nt target sequence) than a user defined threshold is selected as a possible off-target.

Parameters in Command Line. An identifier precedes the value of the parameter, separated by a white space:

  "f1" Precedes Name of the file with the DNA sequences in fasta format to be analysed searching for CRISPR targets.
  "f2" Precedes Name of the file that contains the list of files with the reference genome
  "L0" Precedes Threshold of global mismatches when number of mismatches in seed sequence in zero.
  "L1" Precedes Threshold of global mismatches when number of mismatches in seed sequence in one.
  "ENZ" Precedes identification of endonuclease to be used: Cas9 or Cas12a
  "NAG" Precedes "Yes" or "No", indicating is sequences with PAM "NAG" must be taken into account as possible off-targets 

Using low or high values for L0 and L1 has influence in the determination of the "optimal" candidates. High values increase the number of possible off-targets to be teaken into account, thus all optimal candidates would be very specific and with very low probability of off-target mutations, but it can also happens that no optimal candidates are detected with those parameters.

Examples of command line execution are:

  python ARES_GT_V1.1.py f1 Sequences.txt f2 ChromList.txt L0 5 L1 4 ENZ Cas9 NAG No
 
  python ARES_GT_V1.1.py f1 Sequences.txt f2 ChromList.txt L0 5 L1 4 ENZ Cas12a
 
  python ARES_GT_V1.1.py f1 Genes.txt f2 Contigs.txt L0 4 L1 2 ENZ Cas9 NAG Yes
 
Note that order of parameters has no relevance:

  python ARES_GT_V1.1.py L1 4 f1 Sequences.txt NAG No L0 5 f2 ChromList.txt ENZ Cas9
 
  >>> Examples of input files and output results are available in the ARES-GT folder and in the zip file.

RESULTS

  ARES-GT generates a tabulated text file that is easily exportable to programs as Excel.
  The files contains information about the used files, values of the parameters used in the analysis or total running time.
  A list of all CRISPR candidate targets found in sequences are provided, followed by a second list with all the candidates with only one hit (itself). Then detailed information of each candidate is provided, including all possible off-targets information (as sequences, location and alignment with the candidate CRISPR target). Sense of the sequence are indicated in addition of number of mismatches and the PAM sequence.
  
  As information for all candidates is provided, user can manually identify candidates that have been discarded as optimal (only match itself) due to matching with other genes of a closely related gene family, for instance. This candidates can be very interesting for targeting multiple related genes with only one sgRNA.
  
  Other interesting case is when a very specific location must be targeted. It is possible that with defined parameters (L0 and L1) no optimal candidate had been selected. However, manual evaluation can allow user to select some candidates: if L0 =5 and L1 =4, some possible off-targets can have 4 mismatches distributed in distal sequence (zero mismatches in seed sequence) or 1 mismatch in seed sequence and 2 in distal sequence. Available information about DNA affinity indicates that those cases are very unlikely to be real offtargets "in vivo", moreover if they are not located in the last nucleotides of the target sequence (from PAM).

Example of sequences input file: PDS3_exons.txt
Examples of result files: 2019-01-27_12_29_ALLsgRNAs(PDS3_exons.txt_Cas9).txt
                          2019-01-27_12_30_ALLsgRNAs(PDS3_exons.txt_Cas12a).txt
