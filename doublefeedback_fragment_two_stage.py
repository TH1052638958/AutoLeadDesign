import statistics

import openai
import jsonlines
import time
import pandas as pd
import numpy as np
import random
import json

from rdkit.Chem import Descriptors, QED
from rdkit.Chem import rdMolDescriptors as rdmd
# Set up your API key and model parameters



import openai
import jsonlines
import time
import subprocess
import re
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem import RDConfig
from rdkit.Chem import QED
import os
import sys
# from rdkit.Contrib.SA_Score import sascorer
#change to pre line
# sys.path.append(os.path.join(RDConfig.RDContribDir, 'SA_Score'))
# import sascorer
#change end
from tools.docking import calc_affinity
from tools.tools import get_fragment
from tools.tools import extract_SMILES
import argparse

parser = argparse.ArgumentParser(description='LMLF')
parser.add_argument('--exp-path', type=str, default='6IGX')
parser.add_argument('--protein-name', type=str, default='6igx')
parser.add_argument('--num-molecules', type=int, default=4)
parser.add_argument('--num-generations', type=int, default=100)
parser.add_argument('--threshold', type=float, default=6.5)
parser.add_argument('--temperature', type=float, default=0.6)
parser.add_argument('--model-engine', type=str, default='gpt-4-turbo')
parser.add_argument('--fragment', type=str, default='BRICS')
parser.add_argument('--api-key', type=str, default='sk-JBYywK3F8OYtVWMOC06fB87bE06244B5A5C0726398B6Db64')
#wenwen : https://key.wenwen-ai.com/v1
parser.add_argument('--plantform-url', type=str, default='https://key.wenwen-ai.com/v1')
parser.add_argument('--period', type=int, default=5)
parser.add_argument('--stage-threshold', type=int, default=4)
parser.add_argument('--fragment-process', type=bool, default=True)
parser.add_argument('--cfg-smina', type=str, default='config_smina.yaml')
parser.add_argument('--num-fragment', type=int, default=20)
args = parser.parse_args()

import ruamel.yaml as yaml
from easydict import EasyDict
with open(args.cfg_smina) as f:
    cfg_smina = yaml.load(f, Loader=yaml.Loader)
args.cfg_smina = EasyDict(cfg_smina)


exp_path=args.exp_path

#old docking needed
directory_path = './LMLF/'
protein_name=args.protein_name
num_molecules = args.num_molecules

num_generations = args.num_generations# The number of times to generate new molecules and feed them back into the modelcond
# target_file_path = 'one-box/drd2.jsonl'
target_file_path =os.path.join(exp_path,'init.jsonl')
result_file_path=os.path.join(exp_path,'result.jsonl')
output_file_path = os.path.join(exp_path,'result.csv')
parent_file_path = os.path.join(exp_path,'parent.csv')
from tools.tools import makedir
dir_smina=os.path.join(args.exp_path,'smina')
makedir(dir_smina)
smiles_gen_all=[]
score_docking=[]
par0_list=[]
par1_list=[]
par2_list=[]
#init+gen_good
smiles_all=[]
scores_all=[]

data = []
docking_threshold =args.threshold
unique_molecules = set()
fragment={}
fragment_parent={}
# Generate and feed back new molecules k times
for i in range(1, num_generations):

    
    #Generating molecules
    new_molecules = []
    new_scores=[]

    print("iteration", i)
    smiles_gen_iter=[]
    score_docking_iter=[]

    if (i==1):
        with jsonlines.open(target_file_path) as reader:
            for line in reader:
                if "\n" not in line:
                    # if(line['label']=='1'):
                    new_molecules.append(line['smiles'])
                    unique_molecules.add(line['smiles'])
                    new_scores.append(float(line['score']))


    else:
        # for point in new_target_molecules:
        #     new_molecules.append(point)
        new_molecules=new_target_molecules
        new_scores=new_target_score

    smiles_all+=new_molecules
    scores_all+=new_scores
    for a in range(len(scores_all)):
        if scores_all[a]<0:
            scores_all[a]=0
    for molecule,score in zip(new_molecules,new_scores):
        try:

            fragment_tem=get_fragment(molecule,args.fragment,args.fragment_process)
            if len(fragment_tem)>1:
                for f in fragment_tem:
                    if f in fragment :
                        fragment[f].append(score)
                        fragment_parent[f].append(molecule)
                    else:
                        fragment_score=calc_affinity(f,dir_out=dir_smina,name_protein=protein_name,cfg=args.cfg_smina)
                        fragment[f]=[fragment_score]
                        fragment[f].append(score)
                        fragment_parent[f]=[molecule]
        except Exception as e:
            print('fail to generate fragment')
    df_fragment = pd.DataFrame.from_dict(fragment, orient='index').transpose()
    df_fragment.to_csv(os.path.join(exp_path,'fragment_{}.csv'.format(str(i))), index=False)
    df_fragment_parent = pd.DataFrame.from_dict(fragment_parent, orient='index').transpose()
    df_fragment_parent.to_csv(os.path.join(exp_path,'fragment_parent_{}.csv'.format(str(i))), index=False)

    #fragment filter
    fragment_filted={}
    fragment_list=[]
    fragment_score_list=[]
    score_list_list=[]
    for key,valus in fragment.items():
        fragment_list.append(key)
        fragment_score_list.append(valus[0])
        score_list_list.append(valus[1:])

    if len(fragment_list)<=args.num_fragment:
        indx=range(len(fragment_list))
        print('fragment score:',fragment_score_list)
    else:
        fragment_score_np=np.array(fragment_score_list)
        fragment_np=np.array(fragment_list)
        indx=list(fragment_score_np.argsort()[-args.num_fragment:])
        print('fragment score:', fragment_score_np[indx])
    for idx in indx :
        fragment_filted[fragment_list[idx]]=score_list_list[idx]




    fragment_avg={}
    for key,valus in fragment_filted.items():
        fragment_avg[key]=statistics.mean(valus)
    # fragment_avg_sorted=sorted(fragment_avg.items(),key=lambda d: d[1],reverse=True)
    #fragment
    population=np.array(list(fragment_avg.keys()))
    fragment_avg_score=list(fragment_avg.values())
    for a in range(len(fragment_avg_score)):
        if fragment_avg_score[a]<0:
            fragment_avg_score[a]=0
    weights=np.array(fragment_avg_score)
    #smiles
    population_smiles=np.array(smiles_all)
    weights_smiles=np.array(scores_all)


    qutoe=0
    fragment_concat_set=set()
    num_resampling=0
    while True:
        if num_resampling>=100*args.num_molecules:
            break
        qutoe=qutoe+1
        if i%args.period>=args.stage_threshold:
            sample_smiles = np.random.choice(population_smiles, size=1, replace=False, p=weights_smiles / sum(weights_smiles))
            prompt=f'Generate a novel valid drug-like molecule similar to {sample_smiles[0]} and do not generate any English text'
            sample=[sample_smiles[0],None,None]
        else:
            sample = np.random.choice(population, size=3, replace=False, p=weights / sum(weights))
            fragment_select_tem=sample[0]+sample[1]+sample[2]
            if fragment_select_tem in fragment_concat_set:
                qutoe=qutoe-1

                print('resampling........')
                num_resampling=num_resampling+1
                continue
            else:
                fragment_concat_set.add(fragment_select_tem)


            #print("avg dict:",fragment_avg)
            #print("fragment choised:", sample)
            #change
            prompt = f'Generate a novel valid drug-like molecule SMILES which contains one fragment of [ {sample[0]} , {sample[1]},{sample[2]} ] at least and do not generate any English text.'

        import os
        from openai import OpenAI
        from openai.types.chat import completion_create_params

        os.environ["OPENAI_BASE_URL"] = args.plantform_url
        os.environ["OPENAI_API_KEY"] = args.api_key
        model_engine = args.model_engine  # You can choose a different model if desired
        client = OpenAI()
        time_out=0
        new_mol=None
        while True:
            if new_mol !=None:
                break
            try:
                completion = client.chat.completions.create(
                    model=model_engine,
                    messages=[
                        {"role": "user", "content": prompt}],
                    #n=1,
                    #max_tokens=60,
                    temperature=args.temperature,
                    #stop="!",
                    #user="user"
                    stream=False
                )
                new_mol = json.loads(completion.model_dump_json())["choices"][0]["message"]["content"]
                new_mol=extract_SMILES(new_mol)
            except Exception as e:
                print(f'API ERROR:{e},retrying.......')
                time_out=time_out+1
                if time_out >5:
                    break
                else:
                    continue
        if new_mol==None:
            continue
        if '.' in new_mol:
            continue
        #filte text
        # while True:
        #     filter = r'\'(.*?)\''
        #     import re
        #     new_mol_filted=re.findall(filter,new_mol)
        #     if len(new_mol_filted)!=0:
        #         new_mol=new_mol_filted[0]
        #     else:
        #         break
        try:

            mol = Chem.MolFromSmiles(new_mol)
            if mol is not None :
                #sanitized_mol = Chem.SanitizeMol(mol)
                if new_mol not in unique_molecules:
                    smiles_gen_iter.append(new_mol)
                    unique_molecules.add(new_mol)
                    par0_list.append(sample[0])
                    par1_list.append(sample[1])
                    par2_list.append(sample[2])
                #print("new molecules", new_mol)

        except Exception as e:
            print(f"SMILES Parse Error: {e}. Skipping molecule: {new_mol}")
            continue
        #end
        print(prompt)
        print(new_mol)
        print('num of new mol:',len(smiles_gen_iter))

        if (len(smiles_gen_iter)>=num_molecules) or (qutoe>=2*num_molecules):
            smiles_gen_all+=smiles_gen_iter
            break

    


    for mol in smiles_gen_iter:

        docking_score=calc_affinity(mol,dir_out=dir_smina,name_protein=protein_name,cfg=args.cfg_smina)
        print("docking score", docking_score)
        if(docking_score<0):
            score_docking_iter.append(0)
        else:
            score_docking_iter.append(docking_score)
    score_docking+=score_docking_iter

    
    new_target_molecules=[]
    new_target_score=[]
    for mol, docking_score in zip(smiles_gen_iter, score_docking_iter):
        if (
            docking_score >= docking_threshold

        ):

            new_target_molecules.append(mol)
            new_target_score.append(docking_score)
    print(len(smiles_gen_all),len(score_docking),len(par0_list),len(par1_list),len(par2_list))
    #Put out iterly
    output_iter_path = os.path.join(exp_path, f'result_{i}.csv')
    df = pd.DataFrame({'smiles': smiles_gen_all, 'score': score_docking, 'fragment0': par0_list, 'fragment1': par1_list,
                       'fragment2': par2_list
                       })
    df.to_csv(output_iter_path, index=False)
#Put out finally
df=pd.DataFrame({'smiles':smiles_gen_all,'score':score_docking,'fragment0':par0_list,'fragment1':par1_list,'fragment2':par2_list
                 })

df.to_csv(output_file_path,index=False)

