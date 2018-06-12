#!/usr/bin/env python3
# Name: Andrew Zarzar (azarzar)
# Group Members: None
# 2.1
"""
Supplemental Program. Takes output file from command-motifFinder and makes condensed version in new file. New file
only contains the sequence - percentage within set - number of genes it is found in out of the total genes in the set.
"""
from motifparse import motifparse
import command

"""No calculations, only retrieves information from motifFinder ouput files"""
def main():
    le = 0
    olist = command.fileoutputs[0:3]  # list of file outputs for motifFinder,py
    og = list()  # list of numbers of genes sequence is present in
    for file in olist:
        og.clear()
        motilist = list()
        moti = motifparse(file)
        le = int(moti.returnlength())  # number of gene promoters in file
        for num in moti.returnogtimes():
            og.append(num)  # adds number of genes sequence is present in (from promoters in file)
        for m in moti.getseq():
            motilist.append(m.rstrip())  # adds corresponding sequence

        with open("[1.1]"+file[3:], 'w') as f:
            for i, motif in enumerate(motilist):
                codoncount = {"A": 0, "T": 0, "G": 0, "C": 0}  # used for future GC content implementation
                for char in motif:
                    codoncount[char] += 1
                #p = (0.32**(codoncount.get("A"))) * (0.32**(codoncount.get("T"))) * (0.18**(codoncount.get("G"))) * (0.18**(codoncount.get("C")))
                #n = 500 - len(motif) + 1
                #savg = list()
                average = (int(og[i])) / int(le) * 100
                f.write(str(motif) + " - " + format(average, '.2f') + "% - " + str(og[i]) + "/" + str(le) + "\n")


if __name__ == "__main__":
    main()
