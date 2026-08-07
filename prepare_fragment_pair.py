import numpy as np
from rdkit import Chem
import pandas as pd
import os

import argparse

parser = argparse.ArgumentParser(description='preDiffLinker')
parser.add_argument('--exp-path', type=str, default='/home/th2024/FDB-17/8UOB/')
parser.add_argument('--num', type=int, default=4000)

args = parser.parse_args()

score_file=os.path.join(args.exp_path,'fragment_score.csv')
df=pd.read_csv(score_file)
fragment=df.iloc[:,0]
score=np.array(list(df.iloc[:,1]))
idx=np.array(list(range(len(score))))
weight=score/sum(score)
sample_list=[]
num=0
while True:

    sample_idx = list(np.random.choice(idx, size=2, replace=False, p=weight))
    if (sample_idx in sample_list):
        continue
    else:
        sample_list.append(sample_idx)
        num+=1
    if num>=args.num:
        break
record1=[]
record2=[]
for i in range(len(sample_list)):
    pair=sample_list[i]
    f1=pair[0]
    f2=pair[1]
    record1.append(f1)
    record2.append(f2)
    file1=os.path.join(args.exp_path,'fragment_{}.mol'.format(str(f1)))
    file2 = os.path.join(args.exp_path, 'fragment_{}.mol'.format(str(f2)))
    mol1=Chem.MolFromMolFile(file1)
    mol2 = Chem.MolFromMolFile(file2)
    combined_mol = Chem.CombineMols(mol1, mol2)
    output_file=os.path.join(args.exp_path,'pair_{}.mol'.format(str(i)))
    Chem.MolToMolFile(combined_mol, output_file)
out_record=pd.DataFrame({'f1':record1,'f2':record2})
record_file=os.path.join(args.exp_path,'record.csv')
out_record.to_csv(record_file,index=False)

