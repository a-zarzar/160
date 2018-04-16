#!/usr/bin/env python3 
# Name: Andrew Zarzar (azarzar) 
# Group Members: None

'''
Program parses through the string looking for colons, then seperates info and prints categorized info.
***I know there are easier ways to do the same thing using .split()
Input: fastq string on one line, seperated by colons
Output: fastq information, grouped with respective description/identifier

Assumptions: Colons are used to seperate information. There are no beginning or ending colons.
'''

class FastqString (str):
    ''' Class has 3 methods that are each used to parse through different parts of the string'''
    count = 0
    def parseFirst(self):
        ''' parses through the beginning of the string, looking for seperating semicolons. instantiates and updates count. returns first segment'''
        global count
        count = 0
        data = self
        x = data.index(":", count)
        y = count
        count += x+1
        return self[y:x]
    
    def parseRest(self):
        ''' parses through the middle of the string, looking for seperating semicolons while updating count. returns once a segment is found'''
        global count
        data = self
        x = data.index(":", count)
        y = count
        count = x+1
        return self[y:x]
    
    def parseEnd(self):
        ''' uses count to parse through the end of the string, as data.index wont work for the end of the string'''
        global count
        data = self  
        x = len(data)
        y = count
        count = x+1
        return self[y:x]
    
def main():
    ''' Asks for input, prints the outpus as the input is being parsed for data segments'''
    data = input('Enter FASTQ data: ')
    string = FastqString (data)
    # The methods are used to parse through the string, count is constantly being updated as the methods are used
    #The proper info is being printed
    print ("Instrument = " +string.parseFirst())
    print ("Run ID = "+string.parseRest())
    print ("Flow Cell ID = "+string.parseRest())
    print ("Flow Cell Lane = "+string.parseRest())
    print ("Tile Number = "+string.parseRest())
    print ("X-coord = "+string.parseRest())
    print ("Y-coord = "+string.parseEnd())
    
main()