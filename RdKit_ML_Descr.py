#!/usr/bin/env python


import rdkit
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Draw
from rdkit.Chem.Draw import IPythonConsole
from rdkit.Chem import PandasTools
from rdkit.Chem import Descriptors
from rdkit.ML.Descriptors import MoleculeDescriptors
import pandas as pd
import os
from rdkit import RDConfig



names=['SMILES', 'PDBID']
df=pd.read_csv("Tuttle_pdbseq_smiles", usecols=[0], names=names,sep='\s+',nrows=4000)


# In[4]:


names=['pdb', 'AP']
df2=pd.read_csv('Tuttle_AP_50ns.csv',names=names,nrows=4000,usecols=[1])


frames=[df, df2]
dff = pd.concat(frames, axis=1)
print (dff)


PandasTools.AddMoleculeColumnToFrame(dff, smilesCol='SMILES')


class RDKit_2D:

    def __init__(self, df1):
        self.mols   = [Chem.MolFromSmiles(i) for i in df1['SMILES']]
        self.smiles = [i for i in df1['SMILES']]
        self.ap     = [i for i in df2['AP']]
        #print (self.smiles)
        
    def compute_2Drdkit(self):
        rdkit_2d_desc = []
        calc          = MoleculeDescriptors.MolecularDescriptorCalculator([x[0] for x in Descriptors._descList])
        header        = calc.GetDescriptorNames()
        for i in range(len(self.mols)):
            ds = calc.CalcDescriptors(self.mols[i])
            rdkit_2d_desc.append(ds)
        ML_df = pd.DataFrame(rdkit_2d_desc,   columns=header)
        ML_df.insert(loc=0, column='SMILES',  value=self.smiles)
        ML_df.insert(loc=1, column='AP',      value=self.ap)

        ML_df.to_csv('ML_descriptors_RDKit_2D.csv', index=False)




RDKit_2D(dff).compute_2Drdkit()
