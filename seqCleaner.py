#!/usr/bin/env python3 
# Name: Andrew Zarzar (azarzar) 
# Group Members: None

'''
Reads a DNA string from user input and return a collapsed substring of embedded Ns to: {count}
Any lower case letters are converted to uppercase
Input: a DNA sequence that needs cleaning. (special characters and numbers will not be removed from the final cleaning)
Output: a cleaned DNA sequence that replaces all non-A,T,G,C characters with {numberof non-characters}
*** program will clean any dna string, special characters and numbers will not be removed from the dna string
'''
class DNAstring (str):
    '''returns the length of the string'''
    def length (self):
        return (len(self))
    '''makes dna string all lowercase, only upercases A,T,C,G characters. calls purify method'''
    def sift(self):
        dna = self
        dna = dna.lower()
        dna = dna.replace("a", "A")
        dna = dna.replace("t", "T")
        dna = dna.replace("c", "C")
        dna = dna.replace("g","G")
        newDNA = self.purify(dna) # purify method is called to clean up the dna string
        return newDNA # returns cleaned DNA string
    
    '''creates a list, counts non- ATGC characters and cleans sequence
        Purification removes all non-capitalized letters, replaces it with an updated count enclosed in {}. As the count is updated,
        the old value that was replaced in the previous loop iteration is removed from the list. 
    '''
    def purify(self, sif):
        dna = sif
        l = list(dna) # creats list of dna
        count = 0
        nCount = 1 # the number of non A,T,G,C characters, resets after every A,T,G,C value found
        while count < len(sif): # starts loop that loops through whole input
            c = l[count]
            d = str(l[count-1])
            if c.islower():# if list value at index count isnt A,T,G,C, statement proceeds
                l[count] = "{"+ str(nCount) +"}" #value at index count is turned into "{nCount}"
                if (count!=0) and (d == "{"+ str(nCount-1) +"}"): l[count-1] = "" # statement erases old "{nCount}" value from list
                count = count+1
                nCount = nCount+1
            else:
                nCount = 1 # reset after every A,T,G,C, value found.
                count = count+1      
        return "".join(l)
                
''' Get user DNA data and clean it up.'''      
def main():     
    data = input('Enter DNA data: ')
    thisDNA = DNAstring (data)
    siftDNA = thisDNA.sift() #sifts DNA and calls purify method. Returns the results of the purification
    print ("")
    print("seqCleaner does NOT remove special characters or numbers from the input.\n")
    print ("Cleaned DNA data: " + siftDNA) # prints cleaned DNA
main()



#https://stackoverflow.com/questions/10953189/count-lower-case-characters-in-a-string
    # - provided guidance for implemenation of lowercase character identifier
    
#https://stackoverflow.com/questions/41752946/replacing-a-character-from-a-certain-index
    # - provided guidance for replacing items in a list
    
#https://stackoverflow.com/questions/19970532/how-to-check-a-string-for-a-special-character
    # - provided guidance for implementation of special character identification. (NOT IMPLEMENTED IN CODE ABOVE)
        