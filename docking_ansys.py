import os
import subprocess
from loguru import logger
from rdkit import Chem
from rdkit.Chem import AllChem
import time
from easydict import EasyDict
from tqdm import tqdm
import numpy as np
from tools.docking import calc_affinity
import pandas as pd

if __name__=='__main__':
    # CC1=CC=CN=C1N[C@@H]1CCC2CCCCC21N=Cc1ccc2c(c1)CNC2c1nc(N)c2ccccc2n1
    # smiles = 'C[C@H]1CC[C@H](C2=C3C=CC=CC=C3NC=C2)C[C@H]1[C@H]1CN(C2=C3C=CC=CC=C3NC=C2)C[C@@H]1C'
    '''
    
    smiles_list=['CN1C(=O)C2=CC=CC=C2N=C1NC1=NC2=CC=CC=C2C(N)=N1','O=C(NC1=NC2=CC=CC=C2C(N)=N1)C3CCCN3C(=O)C4=CC=CC=C4',
                 'CC1=CC=C2CNCC2=C1C(=O)NC3=CC=CC(=C3)OC4=NC=C(C5CCNCC5)N=C4','CC1=CC=C2CNCC2=C1C(=O)NC3=CC(=O)C4=C(C=CC=C4O3)C5=CC=CC=C5',
                 'FC1=CN=C2N=C(NC3=NC4=C(C=N3)CCNC4)C(=O)N=C12','ClC1=NC2=C(C=CC=C2)C(=O)N1C1=NC2=CC=CC=C2C(N)=N1',
                 'FC1=CC=C(NC2=CC(I)=C(Cl)C=C2O)C=C1C(=O)NC3CCN(CC3)C4=NC=NC5=C4C=C(Cl)C=C5F','O=C1N(C)C(=NC2=NC3=CC=CC=C3C(N)=N2)N=C2C=CC=CC12',
                 'CN1CCC[C@H]1C2=CC=CC=C2C(=O)NC3=CC=C(CN4CCNCC4)C=C3']
    '''
    smiles_list=['C1=CC=CC=C1C(=O)NCCN1CCN(C2=NC3=CC=CC=C3C(=N2)N)CC1C1=CC=CC=C1']
    # data_df=pd.read_csv('top_8UOB/init.csv')
    # smiles_list=list(data_df.iloc[:,0])
    import ruamel.yaml as yaml
    from easydict import EasyDict
    with open('config_smina_8UOB.yaml') as f:
        cfg_smina = yaml.load(f, Loader=yaml.Loader)
    cfg_smina = EasyDict(cfg_smina)
    for i in range(len(smiles_list)):
        max=0
        max_idx=0
        smiles=smiles_list[i]
        print(f'Docking mol_{i}:{smiles} ................')
        for j in range(10):
            print(f'iter_{j}')
            score=calc_affinity(smiles, dir_out='smina', name_protein='7MJS', cfg=cfg_smina,prefix=f'AAmol_{i}_iter_{j}_',dock_type=0)
            print(score)
            if(score>max):
                max=score
                max_idx=j
            print('max score:',max)
            print('max idx:',max_idx)
        print('-------------------------------------------------------')



