import pandas as pd
import rdkit
from rdkit import Chem
from rdkit.Chem import AllChem
import  os
import argparse
from tools.tools import makedir


parser = argparse.ArgumentParser(description='LMLF')
parser.add_argument('--exp-path', type=str, default='random_8UOB_fragment')

parser.add_argument('--output-dir', type=str, default='init_mol')


args = parser.parse_args()
result_dir=os.path.join(args.exp_path,args.output_dir)
makedir(result_dir)
smile_file=os.path.join(args.exp_path,'init.csv')
data_df=pd.read_csv(smile_file)
smiles=list(data_df.iloc[:,0])
for i in range(len(smiles)):
    sml=smiles[i]
    mol = Chem.MolFromSmiles(sml)
    m2 = Chem.AddHs(mol)
    AllChem.EmbedMolecule(m2)
    m3 = Chem.RemoveHs(m2)
    file_output = os.path.join(result_dir,f'{i}.pdb')
    Chem.MolToPDBFile(m3, file_output)





