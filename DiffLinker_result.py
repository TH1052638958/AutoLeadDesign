import os

import pandas as pd
from rdkit import Chem
exp_path='/home/th2024/DiffLinker/DiffLinker/PRMT5_pair/'
Smiles_list=[]
for i in range(1,2001):
    try:
        file_name=os.path.join(exp_path,'iter_{}'.format(str(i)))
        file_name=os.path.join(file_name,'output_0_pair_{}_.sdf'.format(str(i)))
        mol=Chem.MolFromMolFile(file_name)
        Smiles=Chem.MolToSmiles(mol)
        Smiles_list.append(Smiles)
    except:
        continue
result_file=os.path.join(exp_path,'result.csv')
pd.DataFrame({'smiles':Smiles_list}).to_csv(result_file,index=False)