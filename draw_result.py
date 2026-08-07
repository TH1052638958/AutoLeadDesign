import pandas as pd

from tools.tools import makedir
import argparse
import os
import numpy as np
from rdkit.Chem import Draw
from PIL import Image
from rdkit import Chem







parser = argparse.ArgumentParser(description='LMLF')
parser.add_argument('--result-dir', type=str, default='C_7RPZ_8')
parser.add_argument('--threshold', type=float, default=8)
parser.add_argument('--num-mol', type=int, default=40)
parser.add_argument('--output-dir', type=str, default='C_7RPZ_8')


args = parser.parse_args()
dir_out=os.path.join('result',args.output_dir)
dir_out=os.path.join(dir_out,'target_mol')
makedir(dir_out)
result_file = os.path.join(args.result_dir, 'result.csv')
result_df=pd.read_csv(result_file)
smiles=result_df.iloc[:,0]
scores=result_df.iloc[:,1]
num=0
smiles_filted=[]
scores_filted=[]
idx=[]
for i in range(len(scores)):

    if(scores[i]>=args.threshold):
        num=num+1
        smiles_filted.append(smiles[i])
        scores_filted.append(scores[i])
        idx.append(i)



png_file=os.path.join(dir_out,'result.png')
img = Draw.MolsToGridImage(mols=[Chem.MolFromSmiles(x) for x in smiles_filted],
                                       legends=[str(idx[s])+':'+str(scores_filted[s]) for s in range(len(scores_filted))], molsPerRow=10)
img.save(png_file)

print(num)
