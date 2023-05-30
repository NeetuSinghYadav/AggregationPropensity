#!/usr/bin/env python

import numpy as np


Hydrophobic_aa = {"ALA":"A", "LEU":"L", "VAL":"V", "PHE":"F", "TYR":"Y", "TRP":"W" }
Hydrophilic_aa = {"ASN":"N", "GLN":"Q", "SER":"S", "THR":"T", "HIS":"H"}
Charged_aa     = {"Arg":"R", "Lys":"K", "GLU":"E", "ASP":"D"}
Other_aa       = {"CYS":"C", "Pro":"P", "GLY":"G"}

'''
# # Alpha and beta Structral Propensity Values (taken from https://doi.org/10.1073/pnas.96.22.12524)
# 
# No Values for Proline were listed in the paper, therefore put as 0.0
# 
# ## Hydropbobicity values are taken from (https://www.nature.com/articles/nsb1096-842)
'''

AlphaHelix_Propensity= {"G": 1.24, "A":-0.04, "V":-0.06, "I":-0.26, "L":-0.38,
                        "F":-0.01, "M":-0.09, "W":0.21,  "C":0.57,  "S":0.15,
                        "T": 0.39, "N":0.25,  "Q":-0.02, "H":-0.11, "Y":0.05,
                        "D": 0.27, "E":-0.33, "K":-0.18, "R":-0.30, "P":0.0}


BetaSheet_Propensity=  {"G": 0.76, "A":-0.12, "V":-0.70, "I":-0.77, "L":0.15,
                        "F":-0.67, "M":-0.71, "W":-0.14, "C":-0.63, "S":1.45,
                        "T":-0.70, "N":1.05,  "Q":1.67,  "H":1.34,  "Y":-0.49,
                        "D": 1.12, "E":0.91,  "K":0.29,  "R":0.34,  "P":0.0}


Hydrophobicity      =  {"A":-0.17, "R":-0.81, "N":-0.42, "D":-123, "C":0.24,
                        "Q":-0.58, "E":-2.02, "G":-0.01, "H":-0.17,"I":0.31,
                        "L":0.56,  "K":-0.99, "M":0.23,  "F":1.13, "P":-0.45,
                        "S":-0.13, "T":-0.14, "W":1.85,  "Y":0.94, "V":-0.07}



aa_seq=str("ALRLFNQRS")





aaseq= [i for i in aa_seq]
#print (aaseq)

#aa_list   = ["AAAA", "DRRE", "NQSS", "GSTC", "PQRTTTT"]

#This class caluclates and returns the number of hydrophobic, hydrophilic, chraged and other types of amino acids 
#in the given sequence

class Get_AA_Type:  
    def __init__(self, seq):
        self.seq = seq
        self.hydrophobic_count = 0
        self.hydrophilic_count = 0
        self.other_count       = 0
        self.charged_count     = 0

        #print (self.seq)
        self.xx="".join([str(i) for i in seq])
    def GetAACount(self):
        
        for i in self.seq:
            for k in Hydrophobic_aa.values():
                if i == k:
                    self.hydrophobic_count +=1
                    
            for m in Hydrophilic_aa.values():
                if i == m:
                    self.hydrophilic_count +=1        
            
            for l in Charged_aa.values():
                if i == l:
                    self.charged_count +=1
                
            for n in Other_aa.values():
                if i == n:
                    self.other_count +=1
        
        return [self.xx,self.hydrophobic_count, self.charged_count, self.hydrophilic_count, self.other_count]



class Get_SS_Propensity():
    def __init__(self,seq):
        self.seq           = seq
        self.I_Alpha_count = 0.0
        self.I_Beta_count  = 0.0
        self.xx            = "".join([str(i) for i in seq])

    def alpha_ss_propensity(self):    
        for i in self.seq:
            for key, value in AlphaHelix_Propensity.items():
                if i in key:
                    self.I_Alpha_count +=value
        #print (i, self.I_Alpha_count)            
        return [round(self.I_Alpha_count,3)]
    
    def beta_ss_propensity(self):    
        for i in self.seq:
            for key, value in BetaSheet_Propensity.items():
                if i in key:
                    self.I_Beta_count +=value
        return [round(self.I_Beta_count,3)]    




class Calc_AP_Hydrophobicity():
    def __init__(self,seq):
        self.seq               = str(seq)
        self.P_hydrophobic     = -1.82 
        self.P_Positive_charge =  1.36 
        self.P_negative_charge =  0.95 
        self.P_hydrophilic     =  float(0.10) 
        self.P_other           = -0.61
     
        ss                     = Get_AA_Type(aaseq).GetAACount()
        self.num_hydrophobic   = ss[1]
        self.num_hydrophilic   = ss[2]
        self.num_charge        = ss[3]
       
    def Calc_Ihydr(self):  
        self.Ihydr = 0.0
        for i in self.seq:
            for k in Hydrophobic_aa.values(): 
                if i == k:
                    for key, value in Hydrophobicity.items():
                        if i in key:
                            self.Ihydr +=value-(self.P_hydrophobic*self.num_hydrophobic)
        return (self.Ihydr)
    
    
    def Calc_Ihydrophilicity(self):
        self.Ihydrophil = 0.0
        for i in self.seq:
            for k in Hydrophilic_aa.values(): 
                if i == k:
                    for key, value in Hydrophobicity.items(): # change this with "hydrophilicity"
                        if i in key:
                            self.Ihydrophil +=value-(self.P_hydrophilic*self.num_hydrophilic)
                            #print (i,self.Ihydrophil)
                else:
                    pass
                    #print ("No hydrophilic residues in the seq")
        return (self.Ihydrophil)
                
    

for x in range(len(aa_list)):
    seq             = aa_list[x]
    aaList          = Get_AA_Type(seq).GetAACount()
    I_alpha         = Get_SS_Propensity(seq).alpha_ss_propensity()
    I_beta          = Get_SS_Propensity(seq).beta_ss_propensity()
    Ihydro          = Calc_AP_Hydrophobicity(seq).Calc_Ihydr()
    Ihydrophilicity = Calc_AP_Hydrophobicity(seq).Calc_Ihydrophilicity()
    print (aaList, I_alpha, I_beta, Ihydro,Ihydrophilicity)





def GetAATypes(aa_seq):
    hydrophobic_count = 0
    hydrophilic_count = 0
    charged_count     = 0
    other_count       = 0
    for i in aa_seq:
        for k in Hydrophobic_aa.values():
            if i == k:
                hydrophobic_count+=1

        for l in Charged_aa.values():
            if i == l:
                charged_count    +=1
                
        for m in Hydrophilic_aa.values():
            if i == m:
                hydrophilic_count+=1
           
        for n in Other_aa.values():
            if i == n:
                other_count      +=1
                
    #print ("seq:",aa_seq,
    #       "hydropbobic aa:",hydrophobic_count, "charged aa:",charged_count,
    #       "hydrophilic aa:",hydrophilic_count, "other aa:", other_count) 
    return [aa_seq,hydrophobic_count,charged_count,hydrophilic_count,other_count]



aaList          = Get_AA_Type(aaseq).GetAACount()
Ihydro          = Calc_AP_Hydrophobicity(aaseq).Calc_Ihydr()
Ihydrophilicity = Calc_AP_Hydrophobicity(aaseq).Calc_Ihydrophilicity()
print (aaList, Ihydro,Ihydrophilicity)







