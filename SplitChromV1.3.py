###############
#  Functions  #
###############

import sys


def SplitChrom(filename):
    F = open(filename, "r")
    i = 0
    R = open("Report.txt", "w")
    Genome = F.readlines()
    Lines = ""
    names = []
    for line in Genome:
        line = line.strip()
        if len(line) == 0:
            continue
        if line[0] == ">":
            R.close()
            if line[1:20] in names:
                i = i + 1
                FILE = line[1:20] + "_" + str(i) + ".txt"
            else:
                names.append(line[1:20])
                FILE = line[1:20] + ".txt"
            R = open(FILE, "w")
            R.write(">" + FILE + "\n")
            Lines = Lines + FILE + "\n"
        else:
            R.write(line)
    R.close()
    R = open("ChromosomeList.txt", "w")
    R.write(Lines)
    R.close()


###############
#  Arguments  #
###############

filename = sys.argv[1]

SplitChrom(filename)
