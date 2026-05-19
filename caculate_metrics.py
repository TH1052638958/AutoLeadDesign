
import pandas as pd

from rdkit.Chem import Descriptors, QED
from  rdkit import Chem

import jsonlines
import os
from tools.docking import calc_affinity
from tools.tools import calculate_radscore,calculate_qed_score,obey_lipinski
from tools.tools import makedir
import argparse

parser = argparse.ArgumentParser(description='LMLF')
parser.add_argument('--exp-path', type=str, default='/home/th2024/lmlf/LMLF-main/TSNE_PRMT5_ours1/')
parser.add_argument('--protein-name', type=str, default='7L1G')
parser.add_argument('--init', type=str, default='result.csv')
parser.add_argument('--result', type=str, default='result_others.csv')
parser.add_argument('--threshold', type=float, default=6.5)
parser.add_argument('--cfg-smina', type=str, default='config_smina_PRMT5.yaml')
args = parser.parse_args()


import ruamel.yaml as yaml
from easydict import EasyDict
with open(args.cfg_smina) as f:
    cfg_smina = yaml.load(f, Loader=yaml.Loader)
args.cfg_smina = EasyDict(cfg_smina)
smina_path=os.path.join(args.exp_path,'smina')
makedir(smina_path)
exp_path=args.exp_path
protein_name=args.protein_name
makedir(exp_path)
data=pd.read_csv(os.path.join(exp_path,args.init) )
smiles=list(data.iloc[:,0])
dcoking_scores=[]
mw_scores=[]
logp_scores=[]
qed_scores=[]
rad_scores = []
docking_threshold =args.threshold
mw_threshold = 700
logp_threshold = 6.0
labels=[]
smiles_filted=[]
new_target_molecules=[]
num_lipinski=[]
for smile in smiles:

    dcoking_scores.append(0)
    smiles_filted.append(smile)
    qed_scores.append(calculate_qed_score(smile))
    mw_scores.append(Descriptors.MolWt(Chem.MolFromSmiles(smile)))
    logp_scores.append(Descriptors.MolLogP(Chem.MolFromSmiles(smile)))
    rad_scores.append(calculate_radscore(smile))
    lipinski=obey_lipinski(smile)
    num_lipinski.append(lipinski)
    # radscore_threshold= np.percentile(rad_scores, 75)
for mol, docking_score, qed_score, radscore ,mw,logp in zip(smiles_filted, dcoking_scores, qed_scores, rad_scores,mw_scores,logp_scores):
    if (
        docking_score >= docking_threshold

        # and radscore >= radscore_threshold
        and mw <= mw_threshold
        and logp <= logp_threshold
    ):
        labels.append('1')
        new_target_molecules.append({'smiles': mol, 'label': '1','score':str(docking_score)})
    else:
        labels.append('0')
        new_target_molecules.append({'smiles': mol, 'label': '0','score':str(docking_score)})
with jsonlines.open(os.path.join(exp_path,'init.jsonl') , mode='a') as writer:
        # writer.write("\\n")
    for molecule in new_target_molecules:
        writer.write(molecule)
        writer.write('\n')
pd.DataFrame({'smile':smiles_filted,'dcoking_scores':dcoking_scores,
              'QED':qed_scores,'Rad':rad_scores,'mw':mw_scores,
              'logp':logp_scores,'lipnikin':num_lipinski,'label':labels}).to_csv(os.path.join(exp_path,args.result) ,index=False)
