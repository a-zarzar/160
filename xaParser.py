#!/usr/bin/env python3
# Name: Andrew Zarzar (azarzar)
# Group Members: None
from Bio import AlignIO
from motifparse import motifparse
import scipy.stats as stats
import os
genepromoterlist = ["mes4_up.genes_promoters.txt", "mes4_down.genes_promoters.txt", "mes4_no.change.genes_promoters.txt"]
xlink = [set(), set(), set()]  # xlinked genes [x-up, x-down, x-no change]
alink = [set(), set(), set()]  # autosomal linked genes [ a-up, a-down, a-no]
xdict = [dict(), dict(), dict()]  # dictionary of sequences and their counts in x-linked[up, down, nochange]
adict = [dict(), dict(), dict()]  #||
# up, down, no

def getgene(head):  # gets gene from header
    x = head.index("|") + 1
    genename = head[x:head.index("|", x + 1)]
    return genename


def getchromosome(head):  # gets chromosome from header
    try:
        x = head.index("|") + 1
        gene = head.index("|", x + 1)
        chromosome = head[gene+1:head.index("|", gene+1)]
    except:
        x = head.index("|") + 1
        gene = head.index("|", x + 1)
        chromosome = head[gene + 1:gene+5]
    return str(chromosome)

"""
 Reads the original promoter files records genes as being x-linked or autosomal, storing them in an appropriate list.
 It then parses through the motifFinder outputs and for each sequence, counts how many times it is found in x-linked
 genes, and how many times in autosomal genes. It then updates a dictionary with the sequence and its count.
"""
def main():
    for i, fname in enumerate(genepromoterlist):  # loops through input files (up, down, nochange)
        alignment = AlignIO.read(fname, "fasta")
        for record in alignment:
            if getchromosome(record.id) == "X":
                xlink[i].add(getgene(record.id))  # if on the X, gene is added to xlist at index i
            else:
                alink[i].add(getgene(record.id))  # if autosomal, gene added to alink
    tracker = 0
    for i, fname in enumerate(genepromoterlist):  # loops through input files (up, down, nochange)
        print(fname)
        motilist = list()
        moti = motifparse("[1]"+fname)
        le = moti.getline()  # gets each line from the motifParse output files
        for j in le:
            count = 0
            count2 = 0
            for val in xlink[i]:  # counts number of x-linked (up/down/nochange) genes that motif is present in
                if "'"+str(val)+"'" in j:
                    count += 1
            for val in alink[i]:  # counts number of autosomal (up/down/nochange) genes that motif is present in
                if "'"+str(val)+"'" in j:
                    count2 += 1
            x = {j[0:j.index("[")-1]: count}  # sequence : sequence count
            y = {j[0:j.index("[") - 1]: count2}  #||
            if count != 0:
                xdict[i].update(x)  # updates xdict with sequence and its count in up or down or nochange depending on index
            if count2 != 0:
                adict[i].update(y)  # updates adict
            tracker += 1
            if tracker%200 == 0: print(format((tracker/1553538)*100, '.3f') + "%")
##########################################################conducts statistical analysis in the same fashion as Test.py
    filename = "Outputs/[2]X-linkedSignificant-UP.txt"
    dirname = os.path.dirname(filename)
    if not os.path.exists(dirname):
        os.makedirs(dirname)  # output directory

    with open(filename, 'w') as f:
        f.write(" -Percent Difference- -(1 - hypergeom CDF)-  ----hypergeom PMF---- -sequence- sample, pop. -sample percent- -population percent- " + "\n")
        sorter = list()
        i = 0
        for key, coun in xdict[0].items():
            x = xdict[1].get(key, 1)
            y = xdict[2].get(key, 1)
            z = stats.hypergeom.sf(len(xlink[0]) - 1, 383, int(x) + int(y) + int(coun), int(coun))
            hpd = stats.hypergeom(383, int(x) + int(y) + int(coun), int(len(xlink[0])))
            p = hpd.pmf(coun)
            c = hpd.cdf(coun)
            if p < 0.05 and ((1-c) < 0.05):
                sorter.append([(coun / len(xlink[0])) - ((x + y + coun) / 383), (1-c, float(p), key, coun, (x + y + coun), (coun / len(xlink[0])), ((x + y + coun) / 383))])
            i += 1
            if len(sorter) % 100 == 0: print("#1 / 6 " + format((i / len(xdict[0])) * 100, '.3f') + "%  - " + str(i) + " / " + str(len(xdict[0])))
        print("writing to file..")
        for ls in sorted(sorter, reverse=True):
            f.write(str(ls) + "\n")
#######################################################################
    filename = "Outputs/[2]X-linkedSignificant-DOWN.txt"
    dirname = os.path.dirname(filename)
    if not os.path.exists(dirname):
        os.makedirs(dirname)

    with open(filename, 'w') as f:
        f.write(" -Percent Difference- -(1 - hypergeom CDF)-  ----hypergeom PMF---- -sequence- sample, pop. -sample percent- -population percent- " + "\n")
        sorter = list()
        i = 0
        for key, coun in xdict[1].items():
            x = xdict[0].get(key, 1)
            y = xdict[2].get(key, 1)
            z = stats.hypergeom.sf(len(xlink[0]) -1, 383, int(x) + int(y) + int(coun), int(coun))
            hpd = stats.hypergeom(383, int(x) + int(y) + int(coun), int(len(xlink[1])))
            p = hpd.pmf(coun)
            c = hpd.cdf(coun)
            if p < 0.05 and ((1-c) < 0.05):
                sorter.append([(coun / len(xlink[1])) - ((x + y + coun) / 383), (1-c, float(p), key, coun, (x + y + coun), (coun / len(xlink[1])), ((x + y + coun) / 383))])
            i += 1
            if len(sorter) % 100 == 0: print("#2 / 6" + format((i / len(xdict[1])) * 100, '.3f') + "%  - " + str(i) + " / " + str(len(xdict[1])))
        print("writing to file..")
        for ls in sorted(sorter, reverse=True):
            f.write(str(ls) + "\n")
#########################################################################
    filename = "Outputs/[2]AUTOSOMAL-linkedSignificant-UP.txt"
    dirname = os.path.dirname(filename)
    if not os.path.exists(dirname):
        os.makedirs(dirname)

    with open(filename, 'w') as f:
        f.write(" -Percent Difference- -(1 - hypergeom CDF)-  ----hypergeom PMF---- -sequence- sample, pop. -sample percent- -population percent- " + "\n")
        sorter = list()
        i = 0
        for key, coun in adict[0].items():
            x = adict[1].get(key, 1)
            y = adict[2].get(key, 1)
            hpd = stats.hypergeom(3034, int(x) + int(y) + int(coun), int(len(alink[0])))
            p = hpd.pmf(coun)
            c = hpd.cdf(coun)
            if p < 0.05 and ((1-c) < 0.05):
                sorter.append([(coun / len(alink[0])) - ((x + y + coun) / 3034), (1-c, float(p), key, coun, (x + y + coun), (coun / len(alink[0])), ((x + y + coun) / 3034))])
            i += 1
            if len(sorter) % 100 == 0: print("#3 / 6 " + format((i / len(adict[0])) * 100, '.3f') + "%  - " + str(i) + " / " + str(len(adict[0])))
        print("writing to file..")
        for ls in sorted(sorter, reverse=True):
            f.write(str(ls) + "\n")
##########################################################################
    filename = "Outputs/[2]AUTOSOMAL-linkedSignificant-DOWN.txt"
    dirname = os.path.dirname(filename)
    if not os.path.exists(dirname):
        os.makedirs(dirname)

    with open(filename, 'w') as f:
        f.write(" -Percent Difference- -(1 - hypergeom CDF)-  ----hypergeom PMF---- -sequence- sample, pop. -sample percent- -population percent- " + "\n")
        sorter = list()
        i = 0
        for key, coun in adict[1].items():
            x = adict[0].get(key, 1)
            y = adict[2].get(key, 1)
            hpd = stats.hypergeom(3034, int(x) + int(y) + int(coun), int(len(alink[1])))
            p = hpd.pmf(coun)
            c = hpd.cdf(coun)
            if p < 0.05 and ((1-c) < 0.05):
                sorter.append([(coun / len(alink[1])) - ((x + y + coun) / 3034), (1-c, float(p), key, coun, (x + y + coun), (coun / len(alink[1])), ((x + y + coun) / 3034))])
            i += 1
            if len(sorter) % 100 == 0: print("#4 / 6 " + format((i / len(adict[1])) * 100, '.3f') + "%  - " + str(i) + " / " + str(len(adict[1])))
        print("writing to file..")
        for ls in sorted(sorter, reverse=True):
            f.write(str(ls) + "\n")
#########################################################################################
    filename = "Outputs/[2]AUTOSOMAL-UP-vs x-linked-UP.txt"
    dirname = os.path.dirname(filename)
    if not os.path.exists(dirname):
        os.makedirs(dirname)

    with open(filename, 'w') as f:
        f.write(" -Percent Difference- -(1 - hypergeom CDF)-  ----hypergeom PMF---- -sequence- sample, pop. -sample percent- -population percent- " + "\n")
        sorter = list()
        i = 0
        for key, coun in adict[0].items():
            x = xdict[0].get(key, 1)
            hpd = stats.hypergeom(608, int(x) + int(coun), int(len(alink[0])))
            p = hpd.pmf(coun)
            c = hpd.cdf(coun)
            if p < 0.05 and ((1-c) < 0.05):
                sorter.append([(coun / len(alink[0])) - ((x + coun) / 608), (1-c, float(p), key, coun, (x + coun), (coun / len(alink[0])), ((x + coun) / 608))])
            i += 1
            if len(sorter) % 100 == 0: print("#5 / 6 " + format((i / len(adict[0])) * 100, '.3f') + "%  - " + str(i) + " / " + str(len(adict[0])))
        print("writing to file..")
        for ls in sorted(sorter, reverse=True):
            f.write(str(ls) + "\n")
##############################################################################
    filename = "Outputs/[2]x-linked-UP-vs autosomal-UP.txt"
    dirname = os.path.dirname(filename)
    if not os.path.exists(dirname):
        os.makedirs(dirname)

    with open(filename, 'w') as f:
        f.write(" -Percent Difference- -(1 - hypergeom CDF)-  ----hypergeom PMF---- -sequence- sample, pop. -sample percent- -population percent- " + "\n")
        sorter = list()
        i = 0
        for key, coun in xdict[0].items():
            x = adict[0].get(key, 1)
            hpd = stats.hypergeom(608, int(x) + int(coun), int(len(xlink[0])))
            p = hpd.pmf(coun)
            c = hpd.cdf(coun)
            if p < 0.05 and ((1-c) < 0.05):
                sorter.append([((coun / len(xlink[0])) - ((x + coun) / 608))*100, (1-c, float(p), key, coun, (x + coun), (coun / len(xlink[0]))*100, ((x + coun) / 608)*100)])
            i += 1
            if len(sorter) % 100 == 0: print("#6 / 6 " + format((i / len(xdict[0])) * 100, '.3f') + "%  - " + str(i) + " / " + str(len(xdict[0])))
        print("writing to file..")
        for ls in sorted(sorter, reverse=True):
            f.write(str(ls) + "\n")


if __name__ == "__main__":
    main()
