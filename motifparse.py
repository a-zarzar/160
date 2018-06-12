#!/usr/bin/env python3
# Name: Andrew Zarzar (azarzar)
# Group Members: None
# 3/4
"""
    Contains a number of methods used for parsing FASTA headers, motifFinder output files, and parser output files
    Based on FastAreader provided by David Bernick.
"""
import sys


class motifparse:

    def __init__(self, fname=''):
        self.fname = fname

    def doOpen(self):
        ''' Handle file opens, allowing STDIN.'''
        if self.fname is '':
            return sys.stdin
        else:
            return open(self.fname)

    def getseq(self):  # generator
        with self.doOpen() as fileH:
            for line in fileH:
                x = line.index("[")
                motif = line[0:x]
                yield motif

    def returnlength(self):
        with self.doOpen() as fileH:
            for line in fileH:
                x = line.index("/")
                y = line.index(" ", x)
                return line[x+1:y]

    def returnogtimes(self):  # generator
        with self.doOpen() as fileH:
            for line in fileH:
                x = line.index("/")
                y = line.index("]", x-5)
                yield line[y+1:x]

    def returnlength2(self):
        with self.doOpen() as fileH:
            for line in fileH:
                x = line.index("/")
                y = line.index("-", x)
                return line[x+1:y-1]

    def returnogtimes2(self):  # gets number of genes a sequence is present in when parsing motifFinder output
        with self.doOpen() as fileH:
            for line in fileH:
                x = line.index("/")
                y = line.index("]", x-6)
                yield line[y+1:x]

    def getseq2(self):  # generator
        with self.doOpen() as fileH:
            for line in fileH:
                x = line.index("[")
                motif = line[0:x-1]
                yield motif

    def getline(self):  # generator
        with self.doOpen() as fileH:
            for line in fileH:
                yield line
