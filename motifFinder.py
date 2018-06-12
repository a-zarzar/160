#!/usr/bin/env python3
# Name: Andrew Zarzar (azarzar)
# Group Members: None
# 2
"""
This  Class takes an input and output file as input, and returns all the subsequences in the input file, and the genes
they are found in.
Input: File of FASTA formated headers with promoter sequences of length 500
OutputL File of unique subsequences and the genes that subsequence is found in.
"""
from Bio import AlignIO


class MotifFinder:

    def __init__(self, fname='', allout=''):
        self.fname = fname  # Input File
        self.allout = allout  # Output File
        self.alignment = AlignIO.read(self.fname, "fasta")  # File reader
        self.motiflist = list()  # List of sets. Sets contain all unique subsequences in each promoter
        self.genelist = list()  # Genelist corresponding to the sets found in motiflist
        self.motifdict = dict()  # dictionary of subsequences and all the genes it is found in

    """
        powerSet is initially run from command and takes each promoter sequence from the file adds every subsequence of
        6-mer to 30-mer to a set, ensuring uniqueness. Each set is added to motiflist, with it's gene at a corresponding
        index in genelist  
        Union is called to create a set of all unique subsequences in the whole file by unionizing the whole motiflist
        A dictionary is made with all keys from unionset.
    """
    def powerSet(self):
        count = 0 # tracking method progress
        for record in self.alignment:  # obtains every header in file
            specilist = set()
            for s in range(0, 500):  # promoter is length 500
                for e in range(s + 6, s+31):  # possibe subsequences of 6-mer to 30-mer
                    sub = str(record.seq[s:e])  # every possibe subsequence
                    if len(sub) < 6:
                        break
                    specilist.add(sub)
            self.motiflist.append(specilist)  # adds set to the motiflist
            self.genelist.append(self.getgene(record.id))  # adds gene to the genelist
            count += 1
            avg = (count / len(self.alignment)) * 100 # for printing method progress
            print("Powerset: " + str(count) + "/" + str(len(self.alignment)) + " " + self.fname + " - " + format(avg, '.2f') + "%")
        unionset = self.union(self.motiflist)  # makes union of all sets in motiflist
        self.motifdict = dict.fromkeys(unionset)  # creates dictionary from the unionset values
        self.finder()  # runs finder
        self.writetofile()  # runs writetofile

    """
        After motif dict is complete, sorts data, and outputs to data file.
    """
    def writetofile(self):
        with open(self.allout, 'w') as f:
            moti = list(self.motifdict.items())  # motifdict --> list
            moti.sort(key=lambda x: len(x[1]), reverse=True)  # sort list by genes present in
            moti.sort(key=lambda x: len(x[0]), reverse=True)  # sort list by sequence length
            for key, val in moti:
                out = str((len(val))) + "/" + str(len(self.alignment))
                average = (len(val)) / len(self.alignment) * 100
                if len(val) > 1:  # only adds sequences present in more than 2 genes to file.
                    f.write(str(key) + " " + str(val) + out + " - " + format(average, '.2f') + "%" + "\n")
    """
        Takes union of sets in a list
    """
    @staticmethod
    def union(setlist=list()):
        count = 0
        unionset = set()
        for x in setlist:
            unionset = unionset.union(x)
            count += 1
            print("Union: " + str(count) + "/" + str(len(setlist)))
        return unionset
    """
        Obtains gene name from FASTA header
    """
    def getgene(self, head):
        x = head.index("|")+1
        genename = head[x:head.index("|", x+1)]
        return genename
    """
        Updates motifdict with each unique sequence, and the genes that sequence occurs in.
    """
    def finder(self):
        count = 0  # for method progress
        for key in self.motifdict:  # for every possible subsquence in file, checks for its occurence in each promoter
            genlist = list()
            for i in range(0, len(self.motiflist)):
                if key in self.motiflist[i]:  # if sequence is in promoter, the gene for that promoter ir recorded.
                    genlist.append(self.genelist[i])
            d = {key: genlist}
            self.motifdict.update(d)  # dictionary is updated with sequence, and genes that sequence occurs in
            count += 1
            avg = (count / len(self.motifdict)) * 100  # for printing method progress
            print("Finder: " + str(count) + "/" + str(len(self.motifdict))+ " - " + format(avg, '.2f') + "%")