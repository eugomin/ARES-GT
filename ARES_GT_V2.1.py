###############################################################
# Version 2.0.1                                               #
# Copyright (c) 2020 Eugenio Gomez Minguet.                   #
#                                                             #
# Author : Eugenio Gómez Minguet                              #
# @mail : mingueteg@gmail.com                                 #
###############################################################


from datetime import datetime

import re

import regex

import sys

################
# Diccionaries #
################

codigo = {"A": 0, "C": 1, "G": 2, "T": 3, 'AA': 0, 'AC': 1, 'GT': 11,
          'AG': 2, 'CC': 5, 'CA': 4, 'CG': 6, 'TT': 15, 'GG': 10, 'GC': 9,
          'AT': 3, 'GA': 8, 'TG': 14, 'TA': 12, 'TC': 13, 'CT': 7, 'ACC': 5,
          'ATG': 14, 'AAG': 2, 'AAA': 0, 'ATC': 13, 'AAC': 1, 'ATA': 12,
          'AGG': 10, 'CCT': 23, 'ACT': 7, 'AGC': 9, 'ACA': 4, 'AGA': 8,
          'CAT': 19, 'AAT': 3, 'ATT': 15, 'CTG': 30, 'CTA': 28,
          'CTC': 29, 'CAC': 17, 'ACG': 6, 'CAA': 16, 'AGT': 11, 'CCA': 20,
          'CCG': 22, 'CCC': 21, 'TAT': 51, 'GGT': 43, 'TGT': 59, 'CGA': 24,
          'CAG': 18, 'CGC': 25, 'GAT': 35, 'CGG': 26, 'CTT': 31, 'TGC': 57,
          'GGG': 42, 'TAG': 50, 'GGA': 40, 'TAA': 48, 'GGC': 41, 'TAC': 49,
          'GAG': 34, 'TCG': 54, 'TTA': 60, 'GAC': 33, 'CGT': 27, 'TTT': 63,
          'TCA': 52, 'GCA': 36, 'GTA': 44, 'GCC': 37, 'GTC': 45, 'GCG': 38,
          'GTG': 46, 'TTC': 61, 'GTT': 47, 'GCT': 39, 'TGA': 56, 'TTG': 62,
          'TCC': 53, 'TGG': 58, 'GAA': 32, 'TCT': 55}

codigoi = {"A": "T", "C": "G", "G": "C", "T": "A"}

nt = ("A", "C", "G", "T")


# Sequence antisense-complement
# Transform a sequence in its complement-antisense sequence.
def seq_c(site):
    site_i = site[-1::-1]
    site_c = ""
    for x in site_i:
        if x not in nt:
            site_c = site_c + x
        else:
            y = codigoi[x]
            site_c = site_c + y
    return site_c


########################
# FUNCTION DEFINITIONS #
########################

def ReadFasta(filename):
    F = open(filename, "r")
    sequences = F.read() + "\n"
    F.close()
    v = ""
    m = []
    r = "ok"
    # Regular expressions
    import re
    p = re.compile(r">(?P<name>.+)[^.](?P<seq>[^>]+)")
    p1 = re.compile(r"\s")
    p2 = re.compile(r"[\\/:*?<>|]")
    iterator = p.finditer(sequences)
    # Name and sequence identification
    n = 1
    for seq in iterator:
        sequence = p1.sub("", seq.group("seq")).upper()
        name = p2.sub("", seq.group("name")).upper().replace(" ", "_")
        ACGT = re.search(r"[^ACGT]+", sequence)
        if ACGT:
            v = v + "Ambiguous characters have been found in sequence \"" + name + "\"\
[Sequence number " + str(n) + " in fasta file].\n"
            n = n + 1
        m.append([name, sequence])
        n = n + 1
    v = v + "\n>>>>>>>>>\n\n" + str(len(m)) + " sequences have been imported.\n"
    # To avoid duplicated names
    namelist = []
    for i in range(len(m)):
        name = m[i][0]
        if name in namelist:
            n = 2
            t = 0
            while t == 0:
                nameX = name + "(" + str(n) + ")"
                if nameX not in namelist:
                    t = 1
                    name = nameX
                    m[i][0] = name
                else:
                    n = n + 1
        namelist.append(name)
    if len(m) == 0:
        v = v + "\nNo sequences have been identified, please verify fasta format of sequences file.\n\n"
        r = "end"
    return m, v, r


def ReadChr(Chr):
    F = open(Chr, "r")
    ChrSeq = ""
    Line = F.readlines()
    for line in Line:
        if line[0] == ">":
            continue
        else:
            ChrSeq = ChrSeq + line.replace("\n", "")
    F.close()
    return ChrSeq


def ReadChrlist(filelist):
    F = open(filelist, "r")
    Chrlist = F.readlines()
    F.close()
    Chromosomes = []
    for i in range(len(Chrlist)):
        Chrlist[i] = Chrlist[i].strip()
        if len(Chrlist[i]) == 0:
            continue
        else:
            Chromosomes.append(Chrlist[i])
    return Chromosomes

def ReadChromFASTA(filename):
    F = open(filename, "r")
    i = 0
    R = open("Report.txt", "w")
    Genome = F.readlines()
    Lines = ""
    Contigs = []
    ContigsNames = []
    FILE = "Index"
    C = ""
    for line in Genome:
        line = line.strip()
        if len(line) == 0:
            continue
        if line[0] == ">":
            Contigs.append((FILE, C))
            name = line[1:20].replace(" ","_")
            if name in ContigsNames:
                i = i + 1
                FILE = name + "_" + str(i)
            else:
                FILE = name
            C = ""
            ContigsNames.append(FILE)
        else:
            C = C + line
    Contigs.append((FILE, C))
    Contigs = dict(Contigs[1:])
    ContigsNames.sort()
    return Contigs,ContigsNames

def Identidad(sg4, sg3, SEQ, SEQsense):
    sg4 = sg4.strip()
    Ident = ""
    if sg3 == "+":
        if SEQsense == "+":
            for i in range(len(SEQ)):
                if SEQ[i] == sg4[i]:
                    Ident = Ident + "*"
                else:
                    Ident = Ident + "-"
        elif SEQsense == "-":
            sg4c = seq_c(sg4)
            for i in range(len(SEQ)):
                if SEQ[i] == sg4c[i]:
                    Ident = Ident + "*"
                else:
                    Ident = Ident + "-"
    elif sg3 == "-":
        if SEQsense == "-":
            for i in range(len(SEQ)):
                if SEQ[i] == sg4[i]:
                    Ident = Ident + "*"
                else:
                    Ident = Ident + "-"
        elif SEQsense == "+":
            sg4c = seq_c(sg4)
            for i in range(len(SEQ)):
                if SEQ[i] == sg4c[i]:
                    Ident = Ident + "*"
                else:
                    Ident = Ident + "-"
    return Ident


def alignment(SEQ, SEQsense, Ident):
    A = ""
    if SEQsense == "+":
        for i in range(len(Ident)):
            if Ident[i] == "*":
                A = A + "*"
            else:
                A = A + SEQ[i]
        A = A[:20] + " " + A[20:]
    elif SEQsense == "-":
        for i in range(len(Ident)):
            if Ident[i] == "*":
                A = A + "*"
            else:
                A = A + SEQ[i]
        A = A[:3] + " " + A[3:]
    return A


def PAM(sense, SEQ):
    if sense == "+":
        SEQ = SEQ[:20] + " " + SEQ[20:]
    elif sense == "-":
        SEQ = SEQ[:3] + " " + SEQ[3:]
    return SEQ


def alignment12a(SEQ, SEQsense, Ident):
    A = ""
    if SEQsense == "+":
        for i in range(len(Ident)):
            if Ident[i] == "*":
                A = A + "*"
            else:
                A = A + SEQ[i]
        A = A[:4] + " " + A[4:]
    elif SEQsense == "-":
        for i in range(len(Ident)):
            if Ident[i] == "*":
                A = A + "*"
            else:
                A = A + SEQ[i]
        A = A[:20] + " " + A[20:]
    return A


def PAM12a(sense, SEQ):
    if sense == "+":
        SEQ = SEQ[:4] + " " + SEQ[4:]
    elif sense == "-":
        SEQ = SEQ[:20] + " " + SEQ[20:]
    return SEQ


def MultipleCas9(candidates):
    Targets = []
    Selected = []
    Repeated = []
    for i in range(len(candidates)):
        if candidates[i][3] == "+":
            Targets.append([i, candidates[i][4][:20]])
        elif candidates[i][3] == "-":
            Targets.append([i, seq_c(candidates[i][4][3:])])
    for x in Targets:
        if x[1] in Repeated:
            continue
        else:
            Repeated.append(x[1])
            Multiple = ["", 0]
            for y in Targets:
                if x[0] == y[0]:
                    continue
                else:
                    if x[1] == y[1]:
                        Multiple[0] = Multiple[0] + " = " + candidates[y[0]][0]
                        Multiple[1] = Multiple[1] + 1
                    else:
                        continue
            if Multiple[1] > 0:
                candidates[x[0]][7] = Multiple
                Selected.append(candidates[x[0]])
    return Selected


def MultipleCas9ALL(candidates):
    Targets = []
    Repeated = []
    for i in range(len(candidates)):
        if candidates[i][3] == "+":
            Targets.append([i, candidates[i][4][:20]])
        elif candidates[i][3] == "-":
            Targets.append([i, seq_c(candidates[i][4][3:])])
    for x in Targets:
        if x[1] in Repeated:
            Multiple = [candidates[x[0]][0], -1]
            candidates[x[0]][7] = Multiple
        else:
            Repeated.append(x[1])
            Multiple = ["", 0]
            for y in Targets:
                if x[0] == y[0]:
                    continue
                else:
                    if x[1] == y[1]:
                        Multiple[0] = Multiple[0] + " = " + candidates[y[0]][0]
                        Multiple[1] = Multiple[1] + 1
                    else:
                        continue
            candidates[x[0]][7] = Multiple
    return candidates


def MultipleCas12a(candidates):
    Targets = []
    Selected = []
    Repeated = []
    for i in range(len(candidates)):
        if candidates[i][3] == "+":
            Targets.append([i, candidates[i][4][4:]])
        elif candidates[i][3] == "-":
            Targets.append([i, seq_c(candidates[i][4][:20])])
    for x in Targets:
        if x[1] in Repeated:
            continue
        else:
            Repeated.append(x[1])
            Multiple = ["", 0]
            for y in Targets:
                if x[0] == y[0]:
                    continue
                else:
                    if x[1] == y[1]:
                        Multiple[0] = Multiple[0] + " = " + candidates[y[0]][0]
                        Multiple[1] = Multiple[1] + 1
                    else:
                        continue
            if Multiple[1] > 0:
                candidates[x[0]][7] = Multiple
                Selected.append(candidates[x[0]])
    return Selected


def MultipleCas12aALL(candidates):
    Targets = []
    Repeated = []
    for i in range(len(candidates)):
        if candidates[i][3] == "+":
            Targets.append([i, candidates[i][4][4:]])
        elif candidates[i][3] == "-":
            Targets.append([i, seq_c(candidates[i][4][:20])])
    for x in Targets:
        if x[1] in Repeated:
            Multiple = [candidates[x[0]][0], -1]
            candidates[x[0]][7] = Multiple
        else:
            Repeated.append(x[1])
            Multiple = ["", 0]
            for y in Targets:
                if x[0] == y[0]:
                    continue
                else:
                    if x[1] == y[1]:
                        Multiple[0] = Multiple[0] + " = " + candidates[y[0]][0]
                        Multiple[1] = Multiple[1] + 1
                    else:
                        continue
            candidates[x[0]][7] = Multiple
    return candidates


#############
# Arguments #
#############

# Arguments given when executing program in command line

filename = sys.argv[sys.argv.index("-f1") + 1]  # name of sequences file

ReferenceGenome = sys.argv[sys.argv.index("-f2") + 1]  # Reference Genome in FASTA file

L0 = int(sys.argv[sys.argv.index("-L0") + 1])
if L0 < 0:
    L0 = 0
elif L0 > 9:
    L0 = 9
L0 = L0 + 1
L0NAG = L0 - 1
if L0NAG < 1:
    L0NAG = 1

L1 = int(sys.argv[sys.argv.index("-L1") + 1])
if L1 < 1:
    L1 = 1
elif L1 > 9:
    L1 = 9
L1 = L1 + 1
L1NAG = L1 - 1
if L1NAG < 2:
    L1NAG = 2

Enzyme = sys.argv[sys.argv.index("-ENZ") + 1]

if Enzyme in ["cas9", "Cas9", "CAS9"]:
    Enzyme = "CAS9"
    NAG = sys.argv[sys.argv.index("-NAG") + 1]
    if NAG in ["Y", "y", "YES", "yes", "Yes"]:
        NAG = "Yes"
    else:
        NAG = "No"
else:
    Enzyme = "CAS12a"

OnlyRep = sys.argv[sys.argv.index("-OR") + 1]
if OnlyRep in ["Yes", "Y", "YES", "yes", "y"]:
    OnlyRep = "Y"
else:
    OnlyRep = "N"

print("OR =" + OnlyRep)
#############
# Sequences #
#############

sequences, v, r = ReadFasta(filename)

if r == "end":
    exit()

########################
# Function Execution 1 #
########################

candidates = []

# Candidates internal data format:
# sg = [seq[0],x+1,x+23,"-",SEQ,"","",0,"",0,0,0]
# [ sg[0]  ,    sg[1]      ,   sg[2]     ,     sg[3]   ,  sg[4] ,            sg[5]           ,               sg[6]        , sg[7],          sg[8]             ,
# [Seq name, start position, end position, orientation, sequence, Seed(0) hits mismatches <L0, Seed(1) hits mismatches <L1, Multiple, Seed(2) hits mismatches >L2,
#
#              sg[9]                 ,          sg[10]                    ,              sg[11]                ]
# Counter Seed(0) hits mismatches <L0, Counter Seed(2) hits mismatches >L2, Counter Seed(1) hits mismatches <L1]

if Enzyme == "CAS9":
    for seq in sequences:
        positions = [(a.start()) for a in list(re.finditer('(?=GG)', seq[1][21:]))]
        for x in positions:
            SEQ = seq[1][x:x + 23]
            if "TTTT" in SEQ[:20]:
                continue
            else:
                candidates.append([seq[0], x + 1, x + 23, "+", SEQ, [], [], [], "", 0, 0, 0])
        positions = [(a.start()) for a in list(re.finditer('(?=CC)', seq[1][:-21]))]
        for x in positions:
            SEQ = seq[1][x:x + 23]
            if "AAAA" in SEQ[3:]:
                continue
            else:
                candidates.append([seq[0], x + 1, x + 23, "-", SEQ, [], [], [], "", 0, 0, 0])

    Zeros = len(str(len(candidates)))
    i = 0
    for x in candidates:
        i = i + 1
        x[0] = x[0].strip() + "_" + str(i).zfill(Zeros)

    print("Total number of targets = " + str(len(candidates)))

    if OnlyRep == "Y":
        candidates = MultipleCas9(candidates)
        print("Number of targets with multiple matches = " + str(len(candidates)))
    elif OnlyRep == "N":
        candidates = MultipleCas9ALL(candidates)
        print("OnlyRep = N.")

    now1 = datetime.now()

    FILE = str(now1)[:16].replace(" ", "_").replace(":", "_") + "_ALLsgRNAs(" + filename + "_Cas9).txt"
    R1 = open(FILE, "w")
    R1.write("Analysis of ALL possible sgRNAs targets in sequence(s) from " + filename + "\n\n")
    R1.write("\t\tStarted at " + str(now1)[:19] + "\n\n")
    R1.write("Conditions:\n==========\n\t\tSeed sequence size = 11 nucleotides\n")
    R1.write("\tL0 (NO mismatch in seed sequence)  -> Possible offtarget if global mismatches are less than " + str(
        L0) + "\n")
    R1.write("\tL1 (1 mismatch in seed sequence)   -> Possible offtarget if global mismatches are less than " + str(
        L1) + "\n")
    R1.write("\tTake into account offtargets with \"NAG\" PAM?  " + NAG + "\n")
    R1.write("\tOR = " + OnlyRep + "\n")
    if OnlyRep == "Y":
        R1.write("\n\nWarning! Analysis of ONLY targets that match several query sequences have been selected.\n\n")
    R1.close()

    ########################
    # Function Execution 2 #
    ########################

    Contigs,ContigsNames = ReadChromFASTA(ReferenceGenome)

    for n in ContigsNames:
        ChrN = n
        ChromSeq = Contigs[n]
        for sg in candidates:
            positions = []
            positionsNAG = []
            positionsC = []
            positionsCNAG = []
            if sg[3] == "+":
                seed = sg[4][9:20]
                motif = "(" + seed + "){s<=1}"
                m = regex.finditer(motif, ChromSeq[9:-3], overlapped=True)
                for y in m:
                    x = y.start() + 9
                    if ChromSeq[x + 12:x + 14] in ["GG"]:
                        positions.append(x)
                    if NAG == "Yes":
                        if ChromSeq[x + 12:x + 14] in ["AG"]:
                            positionsNAG.append(x)
                seedC = seq_c(seed)
                motifC = "(" + seedC + "){s<=1}"
                m = regex.finditer(motifC, ChromSeq[3:-9], overlapped=True)
                for y in m:
                    x = y.start()
                    if ChromSeq[x:x + 2] in ["CC"]:
                        positionsC.append(x + 3)
                    if NAG == "Yes":
                        if ChromSeq[x:x + 2] in ["CT"]:
                            positionsCNAG.append(x + 3)
                if len(positions) > 0:
                    for x in positions:
                        SEQ = ChromSeq[x - 9:x + 14]
                        SEQsense = "+"
                        Ident = Identidad(sg[4], sg[3], SEQ, SEQsense)
                        if Ident[9:20].count("-") == 1:
                            if Ident[0:20].count("-") < L1:
                                sg[6].append([ChrN, x - 8, x + 14, "+", SEQ, Ident, Ident[:20].count("-"), "NGG"])
                                sg[11] = sg[11] + 1
                        elif Ident[9:20].count("-") == 0:
                            if Ident[0:20].count("-") < L0:
                                sg[5].append([ChrN, x - 8, x + 14, "+", SEQ, Ident, Ident[:20].count("-"), "NGG"])
                                sg[9] = sg[9] + 1
                if len(positionsC) > 0:
                    for x in positionsC:
                        SEQ = ChromSeq[x - 3:x + 20]
                        SEQsense = "-"
                        Ident = Identidad(sg[4], sg[3], SEQ, SEQsense)
                        if Ident[3:14].count("-") == 1:
                            if Ident[3:].count("-") < L1:
                                sg[6].append([ChrN, x - 2, x + 20, "-", SEQ, Ident, Ident[3:].count("-"), "NGG"])
                                sg[11] = sg[11] + 1
                        elif Ident[3:14].count("-") == 0:
                            if Ident[3:].count("-") < L0:
                                sg[5].append([ChrN, x - 2, x + 20, "-", SEQ, Ident, Ident[3:].count("-"), "NGG"])
                                sg[9] = sg[9] + 1
                if len(positionsNAG) > 0:
                    for x in positionsNAG:
                        SEQ = ChromSeq[x - 9:x + 14]
                        SEQsense = "+"
                        Ident = Identidad(sg[4], sg[3], SEQ, SEQsense)
                        if Ident[9:20].count("-") == 1:
                            if Ident[0:20].count("-") < L1NAG:
                                sg[6].append([ChrN, x - 8, x + 14, "+", SEQ, Ident, Ident[:20].count("-"), "NAG"])
                                sg[11] = sg[11] + 1
                        elif Ident[9:20].count("-") == 0:
                            if Ident[0:20].count("-") < L0NAG:
                                sg[5].append([ChrN, x - 8, x + 14, "+", SEQ, Ident, Ident[:20].count("-"), "NAG"])
                                sg[9] = sg[9] + 1
                if len(positionsCNAG) > 0:
                    for x in positionsCNAG:
                        SEQ = ChromSeq[x - 3:x + 20]
                        SEQsense = "-"
                        Ident = Identidad(sg[4], sg[3], SEQ, SEQsense)
                        if Ident[3:14].count("-") == 1:
                            if Ident[3:].count("-") < L1NAG:
                                sg[6].append([ChrN, x - 2, x + 20, "-", SEQ, Ident, Ident[3:].count("-"), "NAG"])
                                sg[11] = sg[11] + 1
                        elif Ident[3:14].count("-") == 0:
                            if Ident[3:].count("-") < L0NAG:
                                sg[5].append([ChrN, x - 2, x + 20, "-", SEQ, Ident, Ident[3:].count("-"), "NAG"])
                                sg[9] = sg[9] + 1
            elif sg[3] == "-":
                seed = sg[4][3:14]
                motif = "(" + seed + "){s<=1}"
                m = regex.finditer(motif, ChromSeq[3:-9], overlapped=True)
                for y in m:
                    x = y.start()
                    if ChromSeq[x:x + 2] in ["CC"]:
                        positions.append(x + 3)
                    if NAG == "Yes":
                        if ChromSeq[x:x + 2] in ["CT"]:
                            positionsNAG.append(x + 3)
                seedC = seq_c(seed)
                motifC = "(" + seedC + "){s<=1}"
                m = regex.finditer(motifC, ChromSeq[9:-3], overlapped=True)
                for y in m:
                    x = y.start() + 9
                    if ChromSeq[x + 12:x + 14] in ["GG"]:
                        positionsC.append(x)
                    if ChromSeq[x + 12:x + 14] in ["AG"]:
                        if NAG == "Yes":
                            positionsCNAG.append(x)
                if len(positions) > 0:
                    for x in positions:
                        SEQ = ChromSeq[x - 3:x + 20]
                        SEQsense = "-"
                        Ident = Identidad(sg[4], sg[3], SEQ, SEQsense)
                        if Ident[3:14].count("-") == 1:
                            if Ident[3:].count("-") < L1:
                                sg[6].append([ChrN, x - 2, x + 20, "-", SEQ, Ident, Ident[3:].count("-"), "NGG"])
                                sg[11] = sg[11] + 1
                        elif Ident[3:14].count("-") == 0:
                            if Ident[3:].count("-") < L0:
                                sg[5].append([ChrN, x - 2, x + 20, "-", SEQ, Ident, Ident[3:].count("-"), "NGG"])
                                sg[9] = sg[9] + 1
                if len(positionsC) > 0:
                    for x in positionsC:
                        SEQ = ChromSeq[x - 9:x + 14]
                        SEQsense = "+"
                        Ident = Identidad(sg[4], sg[3], SEQ, SEQsense)
                        if Ident[9:20].count("-") == 1:
                            if Ident[0:20].count("-") < L1:
                                sg[6].append([ChrN, x - 8, x + 14, "+", SEQ, Ident, Ident[:20].count("-"), "NGG"])
                                sg[11] = sg[11] + 1
                        elif Ident[9:20].count("-") == 0:
                            if Ident[0:20].count("-") < L0:
                                sg[5].append([ChrN, x - 8, x + 14, "+", SEQ, Ident, Ident[:20].count("-"), "NGG"])
                                sg[9] = sg[9] + 1
                if len(positionsNAG) > 0:
                    for x in positionsNAG:
                        SEQ = ChromSeq[x - 3:x + 20]
                        SEQsense = "-"
                        Ident = Identidad(sg[4], sg[3], SEQ, SEQsense)
                        if Ident[3:14].count("-") == 1:
                            if Ident[3:].count("-") < L1NAG:
                                sg[6].append([ChrN, x - 2, x + 20, "-", SEQ, Ident, Ident[3:].count("-"), "NAG"])
                                sg[11] = sg[11] + 1
                        elif Ident[3:14].count("-") == 0:
                            if Ident[3:].count("-") < L0NAG:
                                sg[5].append([ChrN, x - 2, x + 20, "-", SEQ, Ident, Ident[3:].count("-"), "NAG"])
                                sg[9] = sg[9] + 1
                if len(positionsCNAG) > 0:
                    for x in positionsCNAG:
                        SEQ = ChromSeq[x - 9:x + 14]
                        SEQsense = "+"
                        Ident = Identidad(sg[4], sg[3], SEQ, SEQsense)
                        if Ident[9:20].count("-") == 1:
                            if Ident[0:20].count("-") < L1NAG:
                                sg[6].append([ChrN, x - 8, x + 14, "+", SEQ, Ident, Ident[:20].count("-"), "NAG"])
                                sg[11] = sg[11] + 1
                        elif Ident[9:20].count("-") == 0:
                            if Ident[0:20].count("-") < L0NAG:
                                sg[5].append([ChrN, x - 8, x + 14, "+", SEQ, Ident, Ident[:20].count("-"), "NAG"])
                                sg[9] = sg[9] + 1
            else:
                continue

    now2 = datetime.now()

    R2 = open(FILE, "a")
    if OnlyRep == "N":
        R2.write("\nALL possible Targets found in sequences:\n\n")
    elif OnlyRep == "Y":
        R2.write("\nAnalysed Targets (matching several query sequences):\n\n")

    R2.write("SeqName\tPosStart\tPosEnd\tSense\tSequence(with PAM)\n")
    for seq in candidates:
        hits = seq[9] + seq[11]
        R2.write(seq[0] + "\t" + str(seq[1]) + "\t" + str(seq[2]) + "\t" + seq[3] + "\t" + PAM(seq[3], seq[
            4]) + "\tNumber of Possible OffTargets = " + str(hits - 1) + "\n")

    R2.write("\n\tAnalysis ended at " + str(now1)[:16].replace(" ", "_") + "\n\n")
    R2.write("\t\t\tTotal running time: " + str(now2 - now1) + ".\n\n")

    h = 0
    HITS = ""
    AllTargets = []
    for sg in candidates:
        hits = sg[9] + sg[11]
        if sg[3] == "+":
            AllTargets.append([sg[0], sg[4][:20], seq_c(sg[4][:20]), hits, sg[4][20:], sg[3], sg[7]])
        if sg[3] == "-":
            AllTargets.append([sg[0], sg[4][3:], seq_c(sg[4][3:]), hits, sg[4][:3], sg[3], sg[7]])
        if hits == 1:
            if sg[3] == "+":
                HITS = HITS + sg[0] + "\t" + str(sg[1]) + "\t" + str(sg[2]) + "\t" + sg[3] + "\t" + sg[4][:20] + " " + \
                       sg[4][20:] + "\n"
            elif sg[3] == "-":
                HITS = HITS + sg[0] + "\t" + str(sg[1]) + "\t" + str(sg[2]) + "\t" + sg[3] + "\t" + sg[4][:3] + " " + \
                       sg[4][3:] + "\n"
            h = h + 1

    if OnlyRep == "Y":
        R2.write("\t")
    elif OnlyRep == "N":
        R2.write("Selected Target sequences with Unique Hits:\n\n")
        if h == 0:
            R2.write("SORRY, with this conditions no candidate targets have been found with Unique hits.\n\n")
        else:
            R2.write("SeqName\tPosStart\tPosEnd\tSense\tSequence(with PAM)\n")
            R2.write(HITS)
            HITS = ""

    Multiple = []
    for x in AllTargets:
        if x[6][1] > 0:
            Multiple.append([x[0], x[6][0], x[1], x[3], x[4], x[5], x[6][1]])
    if len(Multiple) > 0:
        R2.write("\nDifferent targets have identical sequence: \n")
        for x in Multiple:
            if x[5] == "+":
                R2.write(x[0] + x[1] + "\t" + x[2] + "  " + "NGG" + "\tNumber of Possible OffTargets = " + str(
                    x[3] - x[6] - 1) + "\n")
            if x[5] == "-":
                R2.write(x[0] + x[1] + "\t" + "CCN" + "  " + x[2] + "\tNumber of Possible OffTargets = " + str(
                    x[3] - x[6] - 1) + "\n")

    R2.write("\nSelected Target sequences with ALL hits:\n")
    R2.write("========================================\n\n")

    from operator import itemgetter

    for sg in candidates:
        R2.write("[[[ TARGET " + sg[0] + " ]]]\n")
        R2.write("\tPosStart\tPosEnd\tSense\tSequence(with PAM)\tReverse\n")
        SEQ = PAM(sg[3], sg[4])
        R2.write("\t" + str(sg[1]) + "\t" + str(sg[2]) + "\t" + sg[3] + "\t" + SEQ + "\t(" + seq_c(SEQ) + ")\n")
        if sg[9] > 0:
            R2.write(">>>><<<<\n>>HITS<<  NO mismatch in Seed sequence (11 nt)\n>>>><<<<\n")
            R2.write("\nChrom\tPosStart\tPosEnd\tSense\tSequence(with PAM)\tSeqIdentity\tMismatches\tPAM\n")
            data = sorted(sg[5], key=itemgetter(6, 0, 1))
            for d in data:
                d[5] = alignment(d[4], d[3], d[5])
                R2.write(d[0] + "\t" + str(d[1]) + "\t" + str(d[2]) + "\t" + d[3] + "\t" + PAM(d[3], d[4]) + "\t" + d[
                    5] + "\t" + str(d[6]) + "\t   " + d[7] + "\n")
            R2.write("\n")
        else:
            R2.write(">>>><<<<\n>>HITS<<  NO mismatch in Seed sequence (11 nt)\n>>>><<<<\n")
            R2.write("\n\t\tNo HITS have been found with this condition.\n\n")
        if sg[11] > 0:
            R2.write(">>>><<<<\n>>HITS<<  1 mismatch in Seed sequence (11 nt)\n>>>><<<<\n")
            R2.write("\nChrom\tPosStart\tPosEnd\tSense\tSequence(with PAM)\tSeqIdentity\tMismatches\tPAM\n")
            data = sorted(sg[6], key=itemgetter(6, 0, 1))
            for d in data:
                d[5] = alignment(d[4], d[3], d[5])
                R2.write(d[0] + "\t" + str(d[1]) + "\t" + str(d[2]) + "\t" + d[3] + "\t" + PAM(d[3], d[4]) + "\t" + d[
                    5] + "\t" + str(d[6]) + "\t   " + d[7] + "\n")
            R2.write("\n")
        else:
            R2.write(">>>><<<<\n>>HITS<<  1 mismatch in Seed sequence (11 nt)\n>>>><<<<\n")
            R2.write("\n\t\tNo HITS have been found with this condition.\n\n")

    R2.close()

    R3 = open("LogInfo.txt", "w")
    R3.write(v)
    R3.close()

#############

elif Enzyme == "CAS12a":
    for seq in sequences:
        positions = [(a.start()) for a in list(re.finditer('(?=AAA)', seq[1][21:]))]
        for x in positions:
            SEQ = seq[1][x:x + 24]
            if "AAAA" in SEQ[:20]:
                continue
            else:
                candidates.append([seq[0], x + 1, x + 24, "-", SEQ, [], [], [], "", 0, 0, 0])
        positions = [(a.start()) for a in list(re.finditer('(?=TTT)', seq[1][:-21]))]
        for x in positions:
            SEQ = seq[1][x:x + 24]
            if "TTTT" in SEQ[4:]:
                continue
            else:
                candidates.append([seq[0], x + 1, x + 24, "+", SEQ, [], [], [], "", 0, 0, 0])

    Zeros = len(str(len(candidates)))
    i = 0
    for x in candidates:
        i = i + 1
        x[0] = x[0].strip() + "_" + str(i).zfill(Zeros)

    print("Total number of targets = " + str(len(candidates)))
    if OnlyRep == "Y":
        candidates = MultipleCas12a(candidates)
        print("Number of targets with multiple matches = " + str(len(candidates)))
    else:
        candidates = MultipleCas12aALL(candidates)

    now1 = datetime.now()

    FILE = str(now1)[:16].replace(" ", "_").replace(":", "_") + "_ALLsgRNAs(" + filename + "_Cas12a).txt"
    R1 = open(FILE, "w")
    R1.write("Analysis of ALL possible sgRNAs targets in sequence(s) from " + filename + "\n\n")
    R1.write("\t\tStarted at " + str(now1)[:19] + "\n\n")
    R1.write("Conditions:\n==========\n\t\tSeed sequence size is 8 nt\n")
    R1.write("\tL0 (NO mismatch in seed sequence)  -> Possible offtarget if global mismatches are less than " + str(
        L0) + "\n")
    R1.write("\tL1 (1 mismatch in seed sequence)   -> Possible offtarget if global mismatches are less than " + str(
        L1) + "\n")
    R1.write("\tOR = " + OnlyRep + "\n")
    if OnlyRep == "Y":
        R1.write("\n\nWarning! Analysis of ONLY targets that match several query sequences have been selected.\n\n")
    R1.close()

    ########################
    # Function Execution 2 #
    ########################

    Contigs,ContigsNames = ReadChromFASTA(ReferenceGenome)

    for n in ContigsNames:
        ChrN = n
        ChromSeq = Contigs[n]
        for sg in candidates:
            positions = []
            positionsC = []
            if sg[3] == "+":
                seed = sg[4][4:12]
                motif = "(" + seed + "){s<=1}"
                m = regex.finditer(motif, ChromSeq[4:-12], overlapped=True)
                for y in m:
                    x = y.start()
                    if ChromSeq[x:x + 3] in ["TTT"]:
                        positions.append(x + 4)
                seedC = seq_c(seed)
                motifC = "(" + seedC + "){s<=1}"
                m = regex.finditer(motifC, ChromSeq[12:-4], overlapped=True)
                for y in m:
                    x = y.start()
                    if ChromSeq[x + 21:x + 24] in ["AAA"]:
                        positionsC.append(x)
                if len(positions) > 0:
                    for x in positions:
                        SEQ = ChromSeq[x - 4:x + 20]
                        SEQsense = "+"
                        Ident = Identidad(sg[4], sg[3], SEQ, SEQsense)
                        if Ident[4:12].count("-") == 1:
                            if Ident[4:24].count("-") < L1:
                                sg[6].append([ChrN, x - 4, x + 19, "+", SEQ, Ident, Ident[4:].count("-"), "NTTT"])
                                sg[11] = sg[11] + 1
                        elif Ident[4:12].count("-") == 0:
                            if Ident[4:24].count("-") < L0:
                                sg[5].append([ChrN, x - 4, x + 19, "+", SEQ, Ident, Ident[4:].count("-"), "NTTT"])
                                sg[9] = sg[9] + 1
                if len(positionsC) > 0:
                    for x in positionsC:
                        SEQ = ChromSeq[x:x + 24]
                        SEQsense = "-"
                        Ident = Identidad(sg[4], sg[3], SEQ, SEQsense)
                        if Ident[12:20].count("-") == 1:
                            if Ident[:20].count("-") < L1:
                                sg[6].append([ChrN, x, x + 23, "-", SEQ, Ident, Ident[:20].count("-"), "NTTT"])
                                sg[11] = sg[11] + 1
                        elif Ident[12:20].count("-") == 0:
                            if Ident[:20].count("-") < L0:
                                sg[5].append([ChrN, x, x + 23, "-", SEQ, Ident, Ident[:20].count("-"), "NTTT"])
                                sg[9] = sg[9] + 1
            elif sg[3] == "-":
                seed = sg[4][12:20]
                motif = "(" + seed + "){s<=1}"
                m = regex.finditer(motif, ChromSeq[12:-4], overlapped=True)
                for y in m:
                    x = y.start()
                    if ChromSeq[x + 21:x + 24] in ["AAA"]:
                        positions.append(x)
                seedC = seq_c(seed)
                motifC = "(" + seedC + "){s<=1}"
                m = regex.finditer(motifC, ChromSeq[4:-12], overlapped=True)
                for y in m:
                    x = y.start()
                    if ChromSeq[x:x + 3] in ["TTT"]:
                        positionsC.append(x + 4)
                if len(positions) > 0:
                    for x in positions:
                        SEQ = ChromSeq[x:x + 24]
                        SEQsense = "-"
                        Ident = Identidad(sg[4], sg[3], SEQ, SEQsense)
                        if Ident[12:20].count("-") == 1:
                            if Ident[:20].count("-") < L1:
                                sg[6].append([ChrN, x - 12, x + 11, "-", SEQ, Ident, Ident[:20].count("-"), "NTTT"])
                                sg[11] = sg[11] + 1
                        elif Ident[9:20].count("-") == 0:
                            if Ident[:20].count("-") < L0:
                                sg[5].append([ChrN, x - 12, x + 11, "-", SEQ, Ident, Ident[:20].count("-"), "NTTT"])
                                sg[9] = sg[9] + 1
                if len(positionsC) > 0:
                    for x in positionsC:
                        SEQ = ChromSeq[x - 4:x + 20]
                        SEQsense = "+"
                        Ident = Identidad(sg[4], sg[3], SEQ, SEQsense)
                        if Ident[4:12].count("-") == 1:
                            if Ident[4:24].count("-") < L1:
                                sg[6].append([ChrN, x - 4, x + 20, "+", SEQ, Ident, Ident[4:].count("-"), "NTTT"])
                                sg[11] = sg[11] + 1
                        elif Ident[4:12].count("-") == 0:
                            if Ident[4:24].count("-") < L0:
                                sg[5].append([ChrN, x - 4, x + 20, "+", SEQ, Ident, Ident[4:].count("-"), "NTTT"])
                                sg[9] = sg[9] + 1
            else:
                continue

    now2 = datetime.now()

    R2 = open(FILE, "a")
    if OnlyRep == "N":
        R2.write("\nALL possible Targets found in sequences:\n\n")
    elif OnlyRep == "Y":
        R2.write("\nAnalysed Targets (matching several query sequences):\n\n")

    R2.write("SeqName\tPosStart\tPosEnd\tSense\tSequence(with PAM)\n")
    for seq in candidates:
        hits = seq[9] + seq[11]
        R2.write(seq[0] + "\t" + str(seq[1]) + "\t" + str(seq[2]) + "\t" + seq[3] + "\t" + PAM12a(seq[3], seq[
            4]) + "\tNumber of Possible OffTargets = " + str(hits - 1) + "\n")

    R2.write("\t\tAnalysis ended at " + str(now1)[:16].replace(" ", "_") + "\n\n")
    R2.write("\t\t\tTotal running time: " + str(now2 - now1) + ".\n\n")

    h = 0
    HITS = ""
    AllTargets = []
    for sg in candidates:
        hits = sg[9] + sg[11]
        if sg[3] == "+":
            AllTargets.append([sg[0], sg[4][4:], seq_c(sg[4][4:]), hits, sg[4][:4], sg[3], sg[7]])
        if sg[3] == "-":
            AllTargets.append([sg[0], sg[4][:20], seq_c(sg[4][:20]), hits, sg[4][20:], sg[3], sg[7]])
        if hits == 1:
            if sg[3] == "+":
                HITS = HITS + sg[0] + "\t" + str(sg[1]) + "\t" + str(sg[2]) + "\t" + sg[3] + "\t" + sg[4][:4] + " " + \
                       sg[4][4:] + "\n"
            if sg[3] == "-":
                HITS = HITS + sg[0] + "\t" + str(sg[1]) + "\t" + str(sg[2]) + "\t" + sg[3] + "\t" + sg[4][:20] + " " + \
                       sg[4][20:] + "\n"
            h = h + 1

    if OnlyRep == "Y":
        R2.write("\t")
    elif OnlyRep == "N":
        R2.write("Selected Target sequences with Unique Hits:\n")
        if h == 0:
            R2.write("\nSORRY, with this conditions no candidate targets have been found with Unique hits.\n\n")
        else:
            R2.write("SeqName\tPosStart\tPosEnd\tSense\tSequence(with PAM)\n")
            R2.write(HITS)
            HITS = ""

    Multiple = []
    for x in AllTargets:
        if x[6][1] > 0:
            Multiple.append([x[0], x[6][0], x[1], x[3], x[4], x[5], x[6][1]])
    if len(Multiple) > 0:
        R2.write("\nDifferent targets have identical sequence: \n")
        for x in Multiple:
            if x[5] == "+":
                R2.write(x[0] + x[1] + "\t" + "TTTN" + "  " + x[2] + "\tNumber of Possible OffTargets = " + str(
                    x[3] - x[6] - 1) + "\n")
            if x[5] == "-":
                R2.write(x[0] + x[1] + "\t" + x[2] + "  " + "NAAA" + "\tNumber of Possible OffTargets = " + str(
                    x[3] - x[6] - 1) + "\n")

    R2.write("\nSelected Target sequences with ALL hits:\n")
    R2.write("========================================\n\n")

    from operator import itemgetter

    for sg in candidates:
        R2.write("[[[ TARGET " + sg[0] + " ]]]\n")
        R2.write("\tPosStart\tPosEnd\tSense\tSequence(with PAM)\tReverse\n")
        SEQ = PAM12a(sg[3], sg[4])
        R2.write("\t" + str(sg[1]) + "\t" + str(sg[2]) + "\t" + sg[3] + "\t" + SEQ + "\t(" + seq_c(SEQ) + ")\n")
        if sg[9] > 0:
            R2.write(">>>><<<<\n>>HITS<<  NO mismatch in Seed sequence (8 nt)\n>>>><<<<\n")
            R2.write("\nChrom\tPosStart\tPosEnd\tSense\tSequence(with PAM)\tSeqIdentity\tMismatches\tPAM\n")
            data = sorted(sg[5], key=itemgetter(6, 0, 1))
            for d in data:
                d[5] = alignment12a(d[4], d[3], d[5])
                R2.write(
                    d[0] + "\t" + str(d[1]) + "\t" + str(d[2]) + "\t" + d[3] + "\t" + PAM12a(d[3], d[4]) + "\t" + d[
                        5] + "\t" + str(d[6]) + "\t   " + d[7] + "\n")
            R2.write("\n")
        else:
            R2.write(">>>><<<<\n>>HITS<<  NO mismatch in Seed sequence (8 nt)\n>>>><<<<\n")
            R2.write("\n\t\tNo HITS have been found with this condition.\n\n")
        if sg[11] > 0:
            R2.write(">>>><<<<\n>>HITS<<  1 mismatch in Seed sequence (8 nt)\n>>>><<<<\n")
            R2.write("\nChrom\tPosStart\tPosEnd\tSense\tSequence(with PAM)\tSeqIdentity\tMismatches\tPAM\n")
            data = sorted(sg[6], key=itemgetter(6, 0, 1))
            for d in data:
                d[5] = alignment12a(d[4], d[3], d[5])
                R2.write(
                    d[0] + "\t" + str(d[1]) + "\t" + str(d[2]) + "\t" + d[3] + "\t" + PAM12a(d[3], d[4]) + "\t" + d[
                        5] + "\t" + str(d[6]) + "\t   " + d[7] + "\n")
            R2.write("\n")
        else:
            R2.write(">>>><<<<\n>>HITS<<  1 mismatch in Seed sequence (8 nt)\n>>>><<<<\n")
            R2.write("\n\t\tNo HITS have been found with this condition.\n\n")

    R2.close()

    R3 = open("LogInfo.txt", "w")
    R3.write(v)
    R3.close()
