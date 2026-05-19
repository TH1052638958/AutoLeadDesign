import pandas as pd

from tools.tools import makedir
import argparse
import os
import numpy as np
from rdkit.Chem import Draw
from PIL import Image
from rdkit import Chem







parser = argparse.ArgumentParser(description='LMLF')
parser.add_argument('--result-dir', type=str, default='random_8UOB_fragment')
parser.add_argument('--threshold', type=float, default=8.5)
parser.add_argument('--num-mol', type=int, default=40)
parser.add_argument('--output-dir', type=str, default='fragment_8UOB')


args = parser.parse_args()
dir_out=os.path.join('result',args.output_dir)
dir_out=os.path.join(dir_out,'target_mol')
makedir(dir_out)
result_file = os.path.join(args.result_dir, 'result.csv')
result_df=pd.read_csv(result_file)
smiles=result_df.iloc[:,0]
scores=result_df.iloc[:,1]
num=0
iter_smiles=[]
iter_scores=[]
iter_idx=[]
for i in range(len(scores)):
    iter=(i//args.num_mol)+1
    iter_dir=os.path.join(dir_out,f'iter_{iter}')
    if(scores[i]>=args.threshold):
        num=num+1
        iter_smiles.append(smiles[i])
        iter_scores.append(scores[i])
        iter_idx.append(i)
    if(i%args.num_mol==args.num_mol-1) and len(iter_smiles)>=1:
        makedir(iter_dir)
        print(iter_idx)
        name=':'.join(str(item) for item in iter_idx)

        png_file=os.path.join(iter_dir,f'{name}.png')
        img = Draw.MolsToGridImage(mols=[Chem.MolFromSmiles(x) for x in iter_smiles],
                                       legends=['%.1f'%s for s in iter_scores], molsPerRow=len(iter_smiles))
        img.save(png_file)
        iter_smiles = []
        iter_scores = []
        iter_idx = []
print(num)




