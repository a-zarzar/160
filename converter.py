#!/usr/bin/env python3 
# Name: Andrew Zarzar (azarzar) 
# Group Members: None

'''
Program will give user instuctions and ask for input. 
    If DNA codons are inputed, the corresponding amino acids are outputed
    If RNA codons are inputed, the corresponding amino acids are outputed
    If a 1 letter amino acid code is inputed, the corresponding 3 letter amino acid code is outputed. Vice-versa
    Program is not case sensitive, and can accept unseperated codons, or codons seperated by spaces or commas. 
    Codons cannot switch between cases
    
    Program differentiates the inputs using a series of specific if-else statements.
    
    *** I took the extra time to design a program that is more robust, and has the ability to translate whole sequences. 
'''
# 3 letter Amino acid to 1 letter Amino acid dictionary
short_AA = {
            'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
            'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N', 
            'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W', 
            'ALA': 'A', 'VAL': 'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'
            }
#creates new dictionary based on short_AA, but swaps items out to allow conversion from 1 AA to 3 AA code.
long_AA = {value:key for key,value in short_AA.items()}

#RNA codon to amino acid dictionary
RNA_codon_table = {
# Second Base
# U             C             A             G
#U
'UUU': 'Phe', 'UCU': 'Ser', 'UAU': 'Tyr', 'UGU': 'Cys',
'UUC': 'Phe', 'UCC': 'Ser', 'UAC': 'Tyr', 'UGC': 'Cys',
'UUA': 'Leu', 'UCA': 'Ser', 'UAA': '---', 'UGA': '---',
'UUG': 'Leu', 'UCG': 'Ser', 'UAG': '---', 'UGG': 'Trp',
#C 
'CUU': 'Leu', 'CCU': 'Pro', 'CAU': 'His', 'CGU': 'Arg',
'CUC': 'Leu', 'CCC': 'Pro', 'CAC': 'His', 'CGC': 'Arg',
'CUA': 'Leu', 'CCA': 'Pro', 'CAA': 'Gln', 'CGA': 'Arg',
'CUG': 'Leu', 'CCG': 'Pro', 'CAG': 'Gln', 'CGG': 'Arg',
#A
'AUU': 'Ile', 'ACU': 'Thr', 'AAU': 'Asn', 'AGU': 'Ser',
'AUC': 'Ile', 'ACC': 'Thr', 'AAC': 'Asn', 'AGC': 'Ser',
'AUA': 'Ile', 'ACA': 'Thr', 'AAA': 'Lys', 'AGA': 'Arg',
'AUG': 'Met', 'ACG': 'Thr', 'AAG': 'Lys', 'AGG': 'Arg',
#G
'GUU': 'Val', 'GCU': 'Ala', 'GAU': 'Asp', 'GGU': 'Gly',
'GUC': 'Val', 'GCC': 'Ala', 'GAC': 'Asp', 'GGC': 'Gly',
'GUA': 'Val', 'GCA': 'Ala', 'GAA': 'Glu', 'GGA': 'Gly',
'GUG': 'Val', 'GCG': 'Ala', 'GAG': 'Glu', 'GGG': 'Gly'
}
# creates DNA codon to amino acid dictionary based on RNA codon table
dnaCodonTable = {key.replace('U','T'):value for key, value in RNA_codon_table.items()}

'''
method is used to loop through multiple rna codons and return amino acids. 
Creates list of every codon from string input.
The list is looped through and replaced with amino acids from dictionary.
The list is joined back into a string
'''
def rnaLoop(str):
    codonList = [ str[i:i+3] for i in range(0, len(str), 3) ]#list of every 3 codons
    translatedList = [RNA_codon_table.get(codonList[i].upper(),"???")+"-" for i in range(0, len(codonList))]#codon--> dictionary
    x = "".join(translatedList)#joins list
    return x[0:len(x)-1]#removes end dash

'''
method is used to loop through multiple dna codons and return amino acids. 
Creates list of every codon from string input.
The list is looped through and replaced with amino acids from dictionary.
The list is joined back into a string
'''
def dnaLoop(str):
    codonList = [ str[i:i+3] for i in range(0, len(str), 3) ]#List of every 3 codons
    translatedList = [dnaCodonTable.get(codonList[i].upper(),"???")+"-" for i in range(0, len(codonList))]#inputs codon into dictionary
    x = "".join(translatedList)#joins list
    return x[0:len(x)-1]#removes end dash

def main():
    print ("Instructions: Converter will convert DNA codons to amino acids, RNA codons to amino acids, and amino acid code to a different amino acid code format  \n" + 
           "Program WILL NOT yield results for a string containing incomplete codons\n"+
           "1. Any string of DNA codons will be converted into an amino acid sequence  \n"+ 
           "2. Any string or RNA codons will be converted into an amino acid sequence \n"+
           "3. Any single amino acid will be converted to its alternate format \n" + "Inputs are NOT case sensitive, and codons can be seperated by spaces,commas, or nothing.\n"+
           "'UNKNOWN' or '???' will be outputted if codon is invalid\n"   )
    data = input("Enter string: ")
    data = data.replace(" ", "") #removes blank spaces from input
    data = data.replace(",","") #removes commas from input
    print("")
    # if input is 3 characters long and is not RNA, then program will search dnaCodon and short_AA dictionaries
    if len(data) ==3 and (data.isupper() or data.islower()) and (('T' in data) or ('t' in data)): print (data.upper()+ " = " + dnaCodonTable.get(data.upper(),short_AA.get(data.upper(),"Unknown")).upper())
    # if input is 3 characters long, then program will search rnaCodon and short_AA dictionaries
    elif len(data)==3 and (data.isupper()or data.islower()): print (data.upper()+ " = " + RNA_codon_table.get(data.upper(),short_AA.get(data.upper(),"Unknown")).upper())
    #satement may be redundant, but will be kept just in case. 
    elif len(data)==3 and not data.isupper() and not data.islower(): print (data.upper()+ " = " + short_AA.get(data.upper(),"Unknown").upper())
    # if input is 1 character long, then program will search long_AA dictionary to get 3 letter amino acids  
    elif len(data)==1: print (data.upper()+ " = " + long_AA.get(data.upper(),"Unknown").upper())
    # if input is a complete set of codons, and is not RNA, then the codons will be converted to 3 letter amino acids
    elif len(data)% 3==0 and (('T' in data) or ('t' in data)): return dnaLoop(data)
    # if input is a complete set of codons, then the codons will be converted to 3 letter amino acids
    elif len(data)% 3==0: return rnaLoop(data)
    # if input isn't a complete set of codons, or a single letter amino acids, then unknown will be returned
    else: print(data+" = UNKNOWN") 
# If dictionary searches turn up negative, then an unknown or ??? value will be returned
main()

# Method rnaLoop based on https://stackoverflow.com/questions/5711452/how-do-i-slice-a-string-every-3-indices
# Method dnaLoop based on rnaLoop