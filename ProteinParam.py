#!/usr/bin/env python3
# Name: Andrew Zarzar (azarzar)
# Group Members: None

class ProteinParam :
    
    '''
    This program is designed to take a protein sequence as an input. The protein sequence must be represented by single-letter 
    amino acids. Many physical-chemical properties of the protein sequence will be calculated. The protein sequence will be 
    read, and the program will print out:
    - The number of amino acids
    - The molecular weight
    - The molar extinction coefficient
    - The mass extinction coefficient
    - The theoretical isoelectic point
    - The amino acid composition
    
    Input: amino acid sequence (in one letter code)
    Output: The various physical-chemical properties outlined above. (Values will be labeled)
    '''
    
# These tables are for calculating:

 #     molecular weight (aa2mw)
    aa2mw = {
        'A': 89.093,  'G': 75.067,  'M': 149.211, 'S': 105.093, 'C': 121.158,
        'H': 155.155, 'N': 132.118, 'T': 119.119, 'D': 133.103, 'I': 131.173,
        'P': 115.131, 'V': 117.146, 'E': 147.129, 'K': 146.188, 'Q': 146.145,
        'W': 204.225,  'F': 165.189, 'L': 131.173, 'R': 174.201, 'Y': 181.189
        }
 # mol. weight of H2O (mwH2O)
    mwH2O = 18.015
 #     absorbance at 280 nm (aa2abs280)
    aa2abs280= {'Y':1490, 'W': 5500, 'C': 125}
 #     pKa of positively charged Amino Acids (aa2chargePos)
    aa2chargePos = {'K': 10.5, 'R':12.4, 'H':6}
 #     pKa of negatively charged Amino acids (aa2chargeNeg)
    aa2chargeNeg = {'D': 3.86, 'E': 4.25, 'C': 8.33, 'Y': 10}
 #     and the constants aaNterm and aaCterm for pKa of the respective terminis
    aaNterm = 9.69
    aaCterm = 2.34

 # becomes a dictionary of the number of each amino acid in the input after aaComposition method runs
    dictionary = {}
 # stringL becomes a string of the cleaned amino acid sequence after aaCount method runs
    stringL = ""
    
    """Provides ProtienParam with the protein attribute. Passes protein to rest of class"""
    def __init__ (self, protein):
        self.protein = protein
        
    """aaCount serves two functions. 1. to count the number of amino acids 2. to create a cleaned amino acid string
    It returns the number of amino acids, but it also cleans a protein list and joins it into a string at the same time.
    """    
    def aaCount (self):
        # list of allowed aa characters to count.
        allowedList = ["A", "C", "D", "E", "F", "G", "H", "I", "L", "K", "M", "N", "P", "Q", "R", "S", "T", "V", "Y", "W"]
        protein = self.protein
        protein = protein.upper() # makes input letters uppercase in case of lower case input.
        l = list(protein) # creates a list of the uppercased protein input.
        count = 0
        aa = 0
        while count < len(protein): # loops through the whole protein sequence.
            c = l[count]
            if c not in allowedList: # if the character at that index in the list isn't an allowed amino acid. It is deleted.
                l[count] = ""
                count += 1
            else:
                aa+= 1 # amino acid count is increased by 1 if character at that index is in the list
                count += 1
      
        global stringL
        stringL =  "".join(l) # joins the cleaned protein list into a string
        return aa # the amino acid count is returned
    
    """pI uses binary search to find the isoelectric point. The precision can be set as a paramater, and affects the precision
    of the pH output to the decimal point. Each search calls the local method _charge_, which returns the net charge. The
    search algorithm tests a new pH in the _charge_ method until a charge of ~ 0 is returned.
    """
    def pI (self, pP): 
        first = 0  # possible pH's from 0 to 14
        last = (14)
        done = False # if true, the while loop below will break. (A suitable pH has been found, or the data has been exhausted)
        precision = pP # Sets precision to have a charge within ?+1 decimal places of zero and a pH accurate to ? decimal places
        p = float("0."+((precision)*"0")+("1")) # translates the users wanted precision to something that can be used
        # if precision is 2, p will equal 0.001, which will result in a pH accurate to 0.01 after the search is completed
        while first <=last and not done:
            mid = (first + last)/2
            x = self._charge_(mid) # tests potential pH for a zero charge using _charge_ method
            if x == 0 or (x<p and x>-p): # if the pH results in a charge sufficiently close to zero (dependent on precision), loop is broken
                done = True
            elif x < 0:
                last = mid - p #if charge is less than 0, then potential pH must be more acidic. Changes test window.
            else: 
                first = mid + p # if charge is greater than 0, then potential pH must be more basic. Changes test window.
        return mid
    
    
    
    """
    **Old, Less efficient pI method. This method tests every possible pH range with a precision of 0.01. Testing potentially
    all 14,000 possible pH values.  
    
    def pI (self):
        first = 0
        last = (14000)
        done = False

        while first <=last and not done:
            x = self._charge_(first/1000)
            if x == 0 or (x<0.01 and x>-0.01):
                done = True
            else: 
                first += 1
   
        return first/1000
     """       
    
    
    
    """aaCompositon tests each value of the allowed amino list against the cleaned protein input string.
    It loops through the list, and counts the instance of each amino acid in the list that are found in the string
    It updates the dictionary outputDict with the amino acid, and the count. Creates a globally available copy of the dictionary
    """
    def aaComposition (self) :
        #list of allowed aa characters
        allowedList = ["A", "C", "D", "E", "F", "G", "H", "I", "L", "K", "M", "N", "P", "Q", "R", "S", "T", "V", "Y", "W"]
        outputDict = {} # new Dictionary. Used to hold the number of each amino acid.
        count = 0
        while count < len(allowedList): #searches the cleaned string from aaCount for the number of every possible amino acid.
            index = allowedList[count] #loops through each character in the allowedList
            x = {index: stringL.count(index)} #sets up new dictionary entry with an amino acid key, and the number of it in the cleaned string
            outputDict.update(x) # updates the dictionary with entry x
            count+=1 # proceeds to the next allowed aa character
        global dictionary
        dictionary = outputDict #creates a copy of outputDict for use outside of the method.
        return (outputDict)
    
    """_charge_ is a private method that is used by the method pI. It returns the net charge of the protein at a specific pH
    and takes a paramater for a pH value. The charge is calculated using an equation that factors the N and C termini, the pH,
    the number of charged and negative amino acids, and their respective pKas.
    """
    def _charge_ (self, p):
        pH = p
        pos = ["K", "R", "H"] # list of positive aa
        neg = ["D", "E", "C", "Y"] # list of negative aa
        count1 = 0
        charge = 0
        while count1 < len(pos): # this loop calculates the positive part of the equation (not including the N terminus)
            index = pos[count1]
            x = 10**ProteinParam.aa2chargePos.get(index) # = 10 to the power of the aa pKa
            y = 10**pH # = 10 to the power of the pH
            z = dictionary.get(index) # the number of the respective aa
            charge +=((x/(x+y))*z) # equation that adds each iteration of the equation for each charged amino acid
            count1 += 1
        x = 10**ProteinParam.aaNterm
        y = 10**pH
        charge += (x/(x+y)) # used to add the N-terimus to the equation seperatly. It would interfere if implemented the same way above.
        
        count2 = 0 
        while count2 < len(neg):# this loop calculates the negative aa part of the equation (not including the C-terminus)
            index = neg[count2]
            x = 10**ProteinParam.aa2chargeNeg.get(index)
            y = 10**pH
            charge -= (y/(x+y))*dictionary.get(index) # equation differs slightly and is instead subtracted from the charge at each iteration
            count2 += 1
        x = 10**ProteinParam.aaCterm
        y = 10**pH
        charge -= (y/(x+y)) # seperate equation for the C-terminus
        return charge
    
    """molarExtinction accepts a cystine paramater. It calculates the extintion coefficient of the protein by using the
    Gill and von Hippel method. It sorts throught the number of tyrosines, tryptophans, and cysteines.
    """
    def molarExtinction (self, cys):
        cystine = cys
        x = 0
        if cystine is not False: cystine = True # under oxidizing conditions
        if cystine is False: x = 1 # if under reducing conditions, cystine doesn't form and cysteine doesn't contribute.
        allowedList = ["Y", "W", "C"] # amino acids required for calculations
        count = 0
        coeff = 0
        while count < (len(allowedList)-x): # finds the data associated with each amino acid, skips cysteine if cystine is False
            index = allowedList[count]
            coeff += (ProteinParam.aa2abs280.get(index)*dictionary.get(index)) # This adds up the number of each aa multiplied by their respective extintion coefficient
            count+=1
        return coeff # returns the protein extintion coefficient
    
    """massExtinction returns the mass extinction coeffectient by using the molar extinction coeffeccient and dividing
    it by the molecular weight of the protein. It also has a cystine paramater. 
    """
    def massExtinction (self,cys):
        cystine = cys
        myMW =  self.molecularWeight() # gets molecular weight
        return self.molarExtinction(cystine) / myMW if myMW else 0.0 # calls for extinction coefficient and divides it by the molecular weight.
    
    """molecularWeight goes through the allowed list of amino acids, and gets the respective value of the molecular weight and the
    number of the amino acid in the input from 2 different dictionaries. It multiplies them together and adds up the weights
    for each amino acid to get the total molecular weight.
    """
    def molecularWeight (self):
        # list of allowed aa characters
        allowedList = ["A", "C", "D", "E", "F", "G", "H", "I", "L", "K", "M", "N", "P", "Q", "R", "S", "T", "V", "Y", "W"]
        count = 0
        weight = 0
        while count < len(allowedList):# loops through each amino acid in the allowed list.
            index = allowedList[count]
            weight += (ProteinParam.aa2mw.get(index)*dictionary.get(index)) # multiplies the number of each aa, and its MW.
            count+=1
        weight = weight - ((len(stringL)-1)*18.01528) # the adjusted weight subtracts the weight of the waters that are released in peptide bond formation
        return weight

# Please do not modify any of the following.  This will produce a standard output that can be parsed

"""The main method creates a ProteinParam object using the input, and prints the useful output by calling the
correct methods. The request for input loops until the program is terminated. """    
import sys
def main():
    inString = input('protein sequence?')
    while inString :
        myParamMaker = ProteinParam(inString)
        myAAnumber = myParamMaker.aaCount()
        myAAcomposition = myParamMaker.aaComposition() # was moved above since some of the methods below rely on its outputs in order to function 
        print ("Number of Amino Acids: {aaNum}".format(aaNum = myAAnumber))
        print ("Molecular Weight: {:.1f}".format(myParamMaker.molecularWeight()))
        print ("molar Extinction coefficient: {:.2f}".format(myParamMaker.molarExtinction(True)))
        print ("mass Extinction coefficient: {:.2f}".format(myParamMaker.massExtinction(True)))
        print ("Theoretical pI: {:.2f}".format(myParamMaker.pI(2)))
        print ("Amino acid composition:")
        keys = list(myAAcomposition.keys())
        keys.sort()
        if myAAnumber == 0 : myAAnumber = 1  # handles the case where no AA are present 
        for key in keys :
            print ("\t{} = {:.2%}".format(key, myAAcomposition[key]/myAAnumber)) # formats the AAcomposition output
        inString = input('protein sequence?') # program will run until terminated

if __name__ == "__main__":
    main()
    
    # Method pI based on http://interactivepython.org/runestone/static/pythonds/SortSearch/TheBinarySearch.html