import argparse
import os.path
import pandas as pd
from tools.tools import get_fragment
from tools.docking import calc_affinity

parser = argparse.ArgumentParser(description='LMLF')
parser.add_argument('--exp-path', type=str, default='fragment_docking')
parser.add_argument('--protein-name',type=str,default='6igx')
parser.add_argument('--fragment-method',type=str,default='BRICS')
parser.add_argument('--fragment-process',type=bool,default=True)
args=parser.parse_args()
data_file=os.path.join(args.exp_path,'init.csv')
data=pd.read_csv(data_file)
smiles=list(data.iloc[:,0])
# smiles=['CCNC(=O)c1ccc(NC[C@@]2(O)CCc3ccccc32)nc1','CN(Cc1nc2ccccc2s1)C(=O)[C@H]1Cc2ccccc2O1','O=C(N[C@H]1CCCNC1=O)c1cc(CCc2ccccc2)ccc1O'
#              ,'Cc1ccc([C@H]2C[C@H]2C(=O)NCc2ccc3nc(O)[nH]c3c2)cc1','O=C(COc1ccc(CO)cc1)N[C@@H]1[C@@H]2Cc3cc(Cl)ccc3[C@@H]21'
#              ,'O=C(CC/C(O)=N/[C@@H]1CCCN=C1O)N1CCSc2ccccc21','NC(=O)C1=NO[C@@H](CNC(=O)N2CCC(c3ccc(O)cc3)CC2)C1']
fragment_all={}
for i in range(len(smiles)):
    fragment = {}
    s=smiles[i]
    score=calc_affinity(s,args.protein_name)
    fragment[s]=[score]
    fragment_tem=list(get_fragment(s,args.fragment_method,args.fragment_process))
    print(fragment_tem)
    for f in fragment_tem:
        score=calc_affinity(f,args.protein_name)
        if f in fragment:
            fragment[f].append(score)
        else:
            fragment[f] = [score]
        if f in fragment_all:
            fragment_all[f].append(score)
        else:
            fragment_all[f] = [score]
    df_fragment = pd.DataFrame.from_dict(fragment, orient='index').transpose()
    df_fragment.to_csv(os.path.join(args.exp_path, f'{i}.csv'), index=False)
df_fragment = pd.DataFrame.from_dict(fragment_all, orient='index').transpose()
df_fragment.to_csv(os.path.join(args.exp_path, 'fragment.csv'), index=False)



