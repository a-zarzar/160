#!/usr/bin/env python3
# Name: Andrew Zarzar (azarzar)
# Group Members: None
# 1
"""
Command Program:

'genepromoterlist' contains the 3 input files
'upmotifs', 'downmotifs', and 'nochangemotifs' are MotifFinder objects that each run a separate input file
Runs:
MotifFinder x3 (Up, Down, No change) --> 3 output files --> Test.py & xaParser.py
parser.py --> 3 output files
Test.py --> 2 output files (Separate Folder)
xaParser.py --> 6 output files (Separate Folder)
"""
from motifFinder import MotifFinder
from threading import Thread
import os
genepromoterlist = ["mes4_up.genes_promoters.txt", "mes4_down.genes_promoters.txt", "mes4_no.change.genes_promoters.txt"]

upmotifs = MotifFinder(genepromoterlist[0], "[1]" + genepromoterlist[0])  # MotifFinder(input file, output file) - UP
downmotifs = MotifFinder(genepromoterlist[1], "[1]" + genepromoterlist[1])  # MotifFinder(input file, output file) - DOWN
nochangemotifs = MotifFinder(genepromoterlist[2], "[1]" + genepromoterlist[2])  # MotifFinder(input file, output file) - No Change

fileoutputs = ["[1]mes4_up.genes_promoters.txt", "[1]mes4_down.genes_promoters.txt", "[1]mes4_no.change.genes_promoters.txt"]

if __name__ == '__main__':
    Thread(target=downmotifs.powerSet()).start()  # runs MotifFinder with downmotifs input
    Thread(target=upmotifs.powerSet()).start()  # runs MotifFinder with upmotifs input
    Thread(target=nochangemotifs.powerSet()).start()  # runs MotifFinder with nochange input
    os.system('parser.py')  # NOT REQUIRED: automatically runs
    os.system('Test.py')  # NOT REQUIRED - UP / DOWN analysis
    os.system('xaParser.py')  # NOT REQUIRED - autosomal / x-linked analysis
