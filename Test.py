#!/usr/bin/env python3
# Name: Andrew Zarzar (azarzar)
# Group Members: None
# 3
"""
Parses the outputs of motfFinder, creates dictionary of each sequence and the time it appears in upregualted genes, down, nochange..
Conducts hypergeometric tests to determine its significance within each sample population
"""
from motifparse import motifparse
import command
import scipy.stats as stats
import os

updict = dict()  # dictionary of upregulated sequences with (# of genes sequence is found in, and total# of genes)
nodict = dict()  # no change
downdict = dict()  # downregulated


def dictmaker(found, length, motifs, file):
    for i, mot in enumerate(motifs):
        d = {mot: (int(found[i]), length)}
        if file == "[1]mes4_up.genes_promoters.txt":
            updict.update(d)
        elif file == "[1]mes4_down.genes_promoters.txt":
            downdict.update(d)
        elif file == "[1]mes4_no.change.genes_promoters.txt":
            nodict.update(d)
"""
Takes each motfFinder output file and gets the number of upregulated genes, and stores the sequence and the number of 
genes it is present in (separate lists). 
Conducts statistical significance test, writes to file
"""
def main():
    olist = command.fileoutputs
    og = list()  # list of number of genes (ints) sequence is present in
    for file in olist:
        print(file)
        og.clear()
        motilist = list()
        moti = motifparse(file)
        le = int(moti.returnlength2())  # number of total genes
        for num in moti.returnogtimes2():  # number of genes sequence present in
            og.append(num)
        for m in moti.getseq2():
            motilist.append(m.rstrip())  # sequence
        dictmaker(og, le, motilist, file)  # sends (number of genes sequence present in, number of genes, sequence list, input file) to dictmaker

    filename = "Outputs/[2]Sig_UP vs down_nochange.txt"
    dirname = os.path.dirname(filename)  # sets output file path
    if not os.path.exists(dirname):      #
        os.makedirs(dirname)             #

    with open(filename, 'w') as f:
        f.write(" -Percent Difference- -(1 - hypergeom CDF)-  ----hypergeom PMF---- -sequence- sample, pop. -sample percent- -population percent- " + "\n")
        sorter = list()
        for key, tup in updict.items():  # key is sequence, tup is number of times sequence appears in upregulated promoters
            x = nodict.get(key, (1, 2651))  # gets number of times sequence appears in null promoters
            y = downdict.get(key, (1, 158))  # gets number of times sequence appears in downregulated promoters
            hpd = stats.hypergeom(3417, int(x[0]) + int(y[0]) + int(tup[0]), int(tup[1]))  # hypergeometric test
            p = hpd.pmf(int(tup[0]))  # probability mass function
            c = hpd.cdf(int(tup[0]))  # cumalitive density function
            if p < 0.05 and ((1-c) < 0.05):  # less than 5% change of data being random
                sorter.append([(tup[0] / tup[1]) - ((x[0] + y[0] + tup[0]) / 3417), (1-c, float(p), key, tup[0], (x[0] + y[0] + tup[0]), (tup[0]/tup[1]), ((x[0] + y[0] + tup[0])/3417))])
        print("writing to file..")
        for ls in sorted(sorter, reverse=True):
            f.write(str(ls) + "\n") # writes to file

    filename = "Outputs/[2]Sig._DOWN vs up_nochange.txt"
    dirname = os.path.dirname(filename)
    if not os.path.exists(dirname):
        os.makedirs(dirname)

    with open(filename, 'w') as f:
        f.write(" -Percent Difference- -(1 - hypergeom CDF)-  ----hypergeom PMF---- -sequence- sample, pop. -sample percent- -population percent- " + "\n")
        sorter2 = list()
        for key, tup in downdict.items():  # times sequence appears in downregulated
            x = nodict.get(key, (1, 2651))  # times appears in null
            y = updict.get(key, (1, 608))  # times appears in upregulated
            hpd = stats.hypergeom(3417, int(x[0]) + int(y[0]) + int(tup[0]), int(tup[1]))  # hypergeometric test
            p = hpd.pmf(int(tup[0]))
            c = hpd.cdf(int(tup[0]))
            if p < 0.05 and ((1-c) < 0.05):
                sorter2.append([(tup[0] / tup[1]) - ((x[0] + y[0] + tup[0]) / 3417), (1 - c, float(p), key, tup[0], (x[0] + y[0] + tup[0]), (tup[0] / tup[1]), ((x[0] + y[0] + tup[0]) / 3417))])
        print("writing to file..")
        for ls in sorted(sorter2, reverse=True):
            f.write(str(ls) + "\n")


if __name__ == "__main__":
    main()
