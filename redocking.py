
import pandas as pd

from rdkit.Chem import Descriptors, QED
from  rdkit import Chem

import jsonlines
import os
from tools.docking import calc_affinity
from tools.tools import calculate_radscore,calculate_qed_score
from tools.tools import makedir
import argparse

parser = argparse.ArgumentParser(description='LMLF')
parser.add_argument('--exp-path', type=str, default='8UOB')
parser.add_argument('--protein-name', type=str, default='8UOB')
parser.add_argument('--file-name', type=str, default='result.csv')
parser.add_argument('--cfg-smina', type=str, default='config_smina.yaml')
args = parser.parse_args()


import ruamel.yaml as yaml
from easydict import EasyDict
with open(args.cfg_smina) as f:
    cfg_smina = yaml.load(f, Loader=yaml.Loader)
args.cfg_smina = EasyDict(cfg_smina)
cfg_smina['autobox_add']='15'

file_path=os.path.join(args.exp_path,args.file_name)
data_df=pd.read_csv(file_path)
smiles=list(data_df.iloc[:,0])
scores=[]
for s in smiles:
    score=calc_affinity(s,dir_out="smina",name_protein=args.protein_name,cfg=args.cfg_smina)
    scores.append(score)
result_df=pd.DataFrame({'smiles':smiles,'score':scores})
result_file=os.path.join(args.exp_path,'result_redocking.csv')
result_df.to_csv(result_file,index=False)
