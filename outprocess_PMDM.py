import os

import pandas as pd
from rdkit import Chem
exp_path='/home/th2024/PMDM/PMDM/data/PRMT5/generate_ref/'
Smiles_list=[]
for i in range(5000):
    file_name=os.path.join(exp_path,'PRMT_{}.sdf'.format(str(i)))
    mol=Chem.MolFromMolFile(file_name)
    Smiles=Chem.MolToSmiles(mol)
    Smiles_list.append(Smiles)
result_file=os.path.join(exp_path,'result.csv')
pd.DataFrame({'smiles':Smiles_list}).to_csv(result_file,index=False)
