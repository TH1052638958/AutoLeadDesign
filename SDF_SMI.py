import os.path

from tools.tools import Sdf2SMILES
import pandas as pd
dir='/home/th2024/PMDM/PMDM/data/8UOB/generate_ref'
out=os.path.join(dir,'result.csv')
SMILES=Sdf2SMILES(dir)
df=pd.DataFrame({'SMILES':SMILES})
df.to_csv(out,index=False)