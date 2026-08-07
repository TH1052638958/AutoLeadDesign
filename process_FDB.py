from rdkit import Chem
import random
import pandas as pd
import os

idx=random.sample(range(1,2000000),2000)

idx_line=0
fragment_list=[]
# 读取 .smi 文件
with open('/home/th2024/FDB-17/FDB-17-fragmentset.smi', 'r') as f:
    for line in f:
        idx_line+=1
        smile = line.strip().split()[0]  # 分割 SMILES 和名称
        if idx_line in idx:
            fragment_list.append(smile)
out=os.path.join('/home/th2024/FDB-17/','sel_2000.csv')
pd.DataFrame({'fragment':fragment_list}).to_csv(out,index=False)