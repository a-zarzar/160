#!/usr/bin/env python3 
# Name: Andrew Zarzar (azarzar) 
# Group Members: None

'''
Read a DNA string from user input and return a collapsed substring of embedded Ns to: {count}

Example: 
 input: AaNNNNNNGTC
output: AA{6}GTC

Any lower case letters are converted to uppercase
'''
class DNAstring (str):
    
    def length (self):
        return (len(self))
    
    def sift(self):
        dna = self
        dna = dna.lower()
        dna = dna.replace("a", "A")
        dna = dna.replace("t", "T")
        dna = dna.replace("c", "C")
        dna = dna.replace("g","G")
        newDNA = self.purify(dna)
        return newDNA
    
    def purify(self, sif):
        dna = sif
        l = list(dna)
        count = 0
        nCount = 1
        while count < len(sif):
            c = l[count]
            d = str(l[count-1])
            if c.islower():
                l[count] = "{"+ str(nCount) +"}"
                if (count!=0) and (d == "{"+ str(nCount-1) +"}"): l[count-1] = ""
                count = count+1
                nCount = nCount+1
            else:
                nCount = 1
                count = count+1      
        return "".join(l)
                
''' Get user DNA data and clean it up.'''      
def main():     
    data = input('Enter DNA data: ')
    thisDNA = DNAstring (data)
    siftDNA = thisDNA.sift()
    print ("")
    print("seqCleaner does NOT remove special characters or numbers from the input.\n")
    print ("Cleaned DNA data: " + siftDNA)
main()



#https://stackoverflow.com/questions/10953189/count-lower-case-characters-in-a-string
    # - provided guidance for implemenation of lowercase character identifier
    
#https://stackoverflow.com/questions/41752946/replacing-a-character-from-a-certain-index
    # - provided guidance for replacing items in a list
    
#https://stackoverflow.com/questions/19970532/how-to-check-a-string-for-a-special-character
    # - provided guidance for inplemenation of special character identification. (NOT IMPLEMENTED IN CODE ABOVE)
        