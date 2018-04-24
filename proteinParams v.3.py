#!/usr/bin/env python3
# Name: Andrew Zarzar (azarzar)
# Group Members: None
# v.3

class ProteinParam :
    
    '''Version 3.0
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
    
    Updates from Version 1.0:
     - Decreased memory requirments by minimizing copies of the input and keeping input in a list.
     - Code for verifying a character is an amino acid was moved from aaCount to aaCOmposition
     - removed redundant and ambiguous code
     
    Updates from Version 2.0:
     - moved creation of composition dictionary to init method.
     - created method attribute self.aaComp of type dict
     - minor implementation changes
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
    
    """Creates the composition dictionary self.aaComp for use by other methods."""
    def __init__ (self, p):
        outputDict = dict.fromkeys(self.aa2mw,0) # creates new Dictionary. Used to hold the number of each amino acid.     
        for count in range(len(p)): #searches the list for the number of every possible amino acid.
            if p[count].upper() in outputDict:
                x = {p[count].upper(): outputDict.get(p[count].upper(),1)+1} #sets up new dictionary entry with an amino acid key, and the number of it in the cleaned string
                outputDict.update(x) # updates the dictionary with entry x
            count+=1 # proceeds to the next allowed aa character
        self.aaComp = outputDict
        
    """aaCount returns the number of amino acids by summing the values of the aa composition dictionary
    """    
    def aaCount (self):
        return sum((self.aaComp).values()) # the amino acid count is returned
    
    """pI uses binary search to find the isoelectric point. The precision can be set as a paramater, and affects the precision
    of the pH output to the decimal point. Each search calls the local method _charge_, which returns the net charge. The
    search algorithm tests a new pH in the _charge_ method until a charge of ~ 0 is returned.
    """
    def pI (self, precision): 
        if self.aaCount() is 0: return 0
        first = 0  # possible pH's from 0 to 14
        last = (14)
        # Sets precision to have a charge within ?+1 decimal places of zero and a pH accurate to ? decimal places
        p = float("0."+((precision)*"0")+("1")) # translates the users wanted precision to something that can be used
        # if precision is 2, p will equal 0.001, which will result in a pH accurate to 0.01 after the search is completed
        while first <=last:
            mid = (first + last)/2
            x = self._charge_(mid) # tests potential pH for a zero charge using _charge_ method
            if x == 0 or (x<p and x>-p): # if the pH results in a charge sufficiently close to zero (dependent on precision), loop is broken
                return mid
            elif x < 0:
                last = mid - p #if charge is less than 0, then potential pH must be more acidic. Changes test window.
            else: 
                first = mid + p # if charge is greater than 0, then potential pH must be more basic. Changes test window.
        return mid
    
    """aaCompositon returns the aaComp dictionary made in the init method.
    """
    def aaComposition (self) :
        return self.aaComp
    
    """_charge_ is a private method that is used by the method pI. It returns the net charge of the protein at a specific pH
    and takes a paramater for a pH value. The charge is calculated using an equation that factors the N and C termini, the pH,
    the number of charged and negative amino acids, and their respective pKas.
    """
    def _charge_ (self, p):
        y = 10**p # = 10 to the power of the pH
        charge = 0
        for aa in self.aa2chargePos.keys(): # this loop calculates the positive part of the equation (not including the N terminus)
            x = 10**self.aa2chargePos.get(aa) # = 10 to the power of the aa pKa
            z = (self.aaComp).get(aa) # the number of the respective aa
            charge +=(x/(x+y))*z # equation that adds each iteration of the equation for each charged amino acid
        x = 10**self.aaNterm
        charge += (x/(x+y)) # used to add the N-terimus to the equation seperatly. It would interfere if implemented the same way above.
        
        for aa in self.aa2chargeNeg.keys():# this loop calculates the negative aa part of the equation (not including the C-terminus)
            x = 10**self.aa2chargeNeg.get(aa)
            z = (self.aaComp).get(aa)
            charge -= (y/(x+y))*z # equation differs slightly and is instead subtracted from the charge at each iteration
        x = 10**self.aaCterm
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
            coeff += (ProteinParam.aa2abs280.get(index)*(self.aaComp).get(index,0)) # This adds up the number of each aa multiplied by their respective extintion coefficient
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
        weight = 0
        for aa in self.aa2mw.keys():# loops through each amino acid in the allowed list.
            weight += (self.aa2mw.get(aa)*self.aaComp.get(aa,0)) # multiplies the number of each aa, and its MW.
        if int(weight) is 0: return weight # handles case when no amino acid is entered
        weight = weight - ((self.aaCount()-1)*18.01528) # the adjusted weight subtracts the weight of the waters that are released in peptide bond formation
        return weight

# Please do not modify any of the following.  This will produce a standard output that can be parsed
    
import sys
def main():
    inString = input('protein sequence?')
    while inString :
        myParamMaker = ProteinParam(inString)
        myAAnumber = myParamMaker.aaCount()
        print ("Number of Amino Acids: {aaNum}".format(aaNum = myAAnumber))
        print ("Molecular Weight: {:.1f}".format(myParamMaker.molecularWeight()))
        print ("molar Extinction coefficient: {:.2f}".format(myParamMaker.molarExtinction(True)))
        print ("mass Extinction coefficient: {:.2f}".format(myParamMaker.massExtinction(True)))
        print ("Theoretical pI: {:.2f}".format(myParamMaker.pI(2)))
        print ("Amino acid composition:")
        myAAcomposition = myParamMaker.aaComposition()
        keys = list(myAAcomposition.keys())
        keys.sort()
        if myAAnumber == 0 : myAAnumber = 1  # handles the case where no AA are present 
        for key in keys :
            print ("\t{} = {:.2%}".format(key, myAAcomposition[key]/myAAnumber))
            
        inString = input('protein sequence?')

if __name__ == "__main__":
    main()