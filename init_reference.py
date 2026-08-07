import argparse
import pandas as pd
import os
from rdkit import Chem
parser = argparse.ArgumentParser(description='LMLF')
parser.add_argument('--out-path', type=str, default='6IGX')
parser.add_argument('--protein-name', type=str, default='6igx')

args = parser.parse_args()
init_file='datasets/init.csv'
df=pd.read_csv(init_file)
smiles=list(df.iloc[:,0])
reference_file_tem='datasets/{}_chainA_ligand'.format(args.protein_name)
if os.path.exists(reference_file_tem+'.pdb'):
    reference_path=reference_file_tem+'.pdb'
    mol=Chem.MolFromPDBFile(reference_path)
else:
    reference_path=reference_file_tem+'.sdf'
    mol=Chem.MolFromMolFile(reference_path)
smi=Chem.MolToSmiles(mol)
smiles.append(smi)
result_file=os.path.join(args.out_path,'init.csv')
pd.DataFrame({'smiles':smiles}).to_csv(result_file,index=False)