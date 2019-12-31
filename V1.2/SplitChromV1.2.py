###############
## Functions ##
###############

def SplitChrom(filename):
    F = open(filename, "r")
    i = 0
    R = open("Report.txt","w")
    Genome = F.readlines()
    Lines = ""
    for line in Genome:
        line = line.strip()
        if line[0] == ">":
            R.close()
            i = i + 1
            FILE = line[1:20]+ "_" + str(i) + ".txt"
            R = open(FILE, "w")
            R.write(">" + FILE + "\n")
            Lines = Lines + FILE + "\n"
        else:
            R.write(line)
    R.close()
    R = open("ChromosomeList.txt","w")
    R.write(Lines)
    R.close()


###############
## Arguments ##
###############

import sys

filename = sys.argv[1]

SplitChrom(filename)


