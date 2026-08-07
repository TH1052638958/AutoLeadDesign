import pandas as pd
import rdkit
from rdkit import Chem
result_file='/home/th2024/lmlf/LMLF-main/TSNE_8UOB_ours/result_19.csv'

df=pd.read_csv(result_file)
smiles=list(df.iloc[:,0])
scores=list(df.iloc[:,1])
f0=list(df.iloc[:,2])
f1=list(df.iloc[:,3])
f2=list(df.iloc[:,4])
print(scores[739])
smiles.insert(740,'CC1=C(NC(=O)C2=CC=CC(C)=C2C2=NC(N)=CC=C2)C=CC=C1N')
scores.insert(740,9.2)
f0.insert(740,'NC1=C(C)C(N)=CC=C1')
f1.insert(740,'C1=NC(N)=CC=C1')
f2.insert(740,'C1=CC(C)=CC=C1')
pd.DataFrame({'smiles': smiles, 'score': scores, 'fragment0': f0, 'fragment1': f1,
                       'fragment2': f2
                       }).to_csv('/home/th2024/lmlf/LMLF-main/TSNE_8UOB_ours/result.csv',index=False)