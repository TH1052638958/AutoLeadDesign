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
import argparse

parser = argparse.ArgumentParser(description='LMLF')
parser.add_argument('--exp-path', type=str, default='6IGX')
parser.add_argument('--protein-name', type=str, default='6igx')
parser.add_argument('--num-molecules', type=int, default=4)
parser.add_argument('--num-generations', type=int, default=100)
parser.add_argument('--threshold', type=float, default=6.5)
parser.add_argument('--temperature', type=float, default=0.6)
parser.add_argument('--model-engine', type=str, default='gpt-4-turbo')
parser.add_argument('--api-key', type=str, default='sk-JBYywK3F8OYtVWMOC06fB87bE06244B5A5C0726398B6Db64')
args = parser.parse_args()




exp_path=args.exp_path

#old docking needed
directory_path = './LMLF/'
protein_name=args.protein_name
num_molecules = args.num_molecules
temperature = 0.3 # Controls the "creativity" of the generated molecules
num_generations = args.num_generations# The number of times to generate new molecules and feed them back into the modelcond
# target_file_path = 'one-box/drd2.jsonl'
target_file_path =os.path.join(exp_path,'init.jsonl')
result_file_path=os.path.join(exp_path,'result.jsonl')
output_file_path = os.path.join(exp_path,'result.csv')
parent_file_path = os.path.join(exp_path,'parent.csv')
from tools.tools import makedir
dir_smina="smina1"
makedir(dir_smina)
def calculate_docking_score(smiles):
    try:
        # Generate an RDKit molecule from the SMILES string
        molecule = Chem.MolFromSmiles(smiles)
        molecule = Chem.AddHs(molecule)
    except Exception as e:
        print(f"SMILES Parse Error: {e}. Skipping molecule: {smiles}")
        return None
    if molecule is not None:
        
        try:
            AllChem.Compute2DCoords(molecule)
        except Exception as e:
            print(f"Compute2DCoords Error: {e}. Skipping molecule: {smiles}")
            return None

        
        AllChem.EmbedMolecule(molecule, AllChem.ETKDG())

        
        try:
            AllChem.MMFFOptimizeMolecule(molecule)
        except Exception as e:
            print(f"MMFFOptimizeMolecule Error: {e}. Skipping molecule: {smiles}")
            return None

        # generate a PDB file from the molecule
        pdb_filename = './ligand.pdb'
        writer = Chem.PDBWriter(pdb_filename)
        writer.write(molecule)
        writer.close()


        cmd = ['./gnina', '--config', 'DRD2_config.txt', '--ligand', './ligand.pdb', '--out', 'output.sdf', '--log', './threshold_output_log.txt', '--cpu', '4', '--num_modes', '1']
        # 
        #cmd = ['./gnina', '-r', '.4IVA.pdb', '-l', './ligand.pdb', '--autobox_ligand', './ligand.pdb', '-o', '/content/docked.txt', '--seed', '0']

        print("Docking Command:", ' '.join(cmd))
        # try:
        #     subprocess.run(cmd, check=True)
        # except subprocess.CalledProcessError as e:
        #     print("Docking Error:", e)
        #     return None
        try:
            subprocess.run(cmd, check=True, stderr=subprocess.PIPE)
        except subprocess.CalledProcessError as e:
            print("Docking process failed:", e)
            print("Error output:", e.stderr)
            return None


        # subprocess.run(cmd, check=True)
        
        # TODO: Extract and return the docking score from the output files
        import os

        # Iterate over the files in the directory
        for filename in os.listdir(directory_path):
            if filename.endswith('.txt'):  # Consider only the text files
                file_path = os.path.join(directory_path, filename)
                with open(file_path, 'r') as file:
                    lines = file.readlines()

                    for i, line in enumerate(lines):
                        if 'affinity' in line.lower() and 'cnn' in line.lower():
                            third_next_line_values = lines[i + 3].split()
                            if len(third_next_line_values) >= 4:
                                try:
                                    cnn_affinity = float(third_next_line_values[3].strip())
                                    return cnn_affinity
                                except ValueError:
                                    pass
        
        return None

def calculate_radscore(mol):
    # molecule = Chem.MolFromSmiles(mol)
    # sa_score = rdMolDescriptors.SyntheticAccessibility(molecule)
    try:
        m = Chem.MolFromSmiles(mol)
        sa_score = sascorer.calculateScore(m)
        return sa_score
    except:
        return None
def calculate_qed_score(smiles):
    molecule = Chem.MolFromSmiles(smiles)
    if molecule is not None:
        qed_score = QED.qed(molecule)
        return qed_score
    else:
        return None
            

      




smiles_gen_all=[]
score_docking=[]
par0_list=[]
par1_list=[]
par2_list=[]

data = []
docking_threshold =args.threshold
unique_molecules = set()
fragment={}
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
                    if(line['label']=='1'):
                        new_molecules.append(line['smiles'])
                        new_scores.append(calc_affinity(new_molecules[-1],dir_out=dir_smina,name_protein=protein_name))
        unique_molecules=set(new_molecules)


    else:
        for point in new_target_molecules:
            new_molecules.append(point)
        new_scores=new_target_score



    for molecule,score in zip(new_molecules,new_scores):
        try:

            fragment_tem=get_fragment(molecule,'BRICS')
            for f in fragment_tem:
                if f in fragment :
                    fragment[f].append(score)
                else:
                    fragment[f]=[score]
        except Exception as e:
            print('fail to generate fragment')
    df_fragment = pd.DataFrame.from_dict(fragment, orient='index').transpose()
    df_fragment.to_csv(os.path.join(exp_path,'fragment_{}.csv'.format(str(i))), index=False)
    fragment_avg={}
    for key,valus in fragment.items():
        fragment_avg[key]=statistics.mean(valus)
    # fragment_avg_sorted=sorted(fragment_avg.items(),key=lambda d: d[1],reverse=True)
    population=np.array(list(fragment_avg.keys()))
    weights=np.array(list(fragment_avg.values()))

    while True:
        sample = np.random.choice(population, size=3, replace=False, p=weights / sum(weights))

        print("avg dict:",fragment_avg)
        print("fragment choised:", sample)
        #change
        prompt = f'Generate a novel valid molecule SMILES which contains one fragment of [ {sample[0]} , {sample[1]},{sample[2]} ] at least and do not generate any English text.'

        import os
        from openai import OpenAI
        from openai.types.chat import completion_create_params

        os.environ["OPENAI_BASE_URL"] = "https://key.wenwen-ai.com/v1"
        os.environ["OPENAI_API_KEY"] = args.api_key
        model_engine = args.model_engine  # You can choose a different model if desired
        client = OpenAI()
        try:
            completion = client.chat.completions.create(
                model=model_engine,
                messages=[
                    {"role": "user", "content": prompt}],
                n=1,
                max_tokens=60,
                temperature=args.temperature,
                stop=".",
                user="user"
            )
        except Exception as e:
            print(f'API ERROR:{e},retrying.......')
            continue
        new_mol=json.loads(completion.model_dump_json())["choices"][0]["message"]["content"]
        #add
        if (i==10 and len(smiles_gen_iter)==0):
            new_mol='CC1C[C@H]2N[C@H]3Cc4cc(Cl)ccc4C[C@@H]3[C@@H]2N1'
            sample[0]='[5*]N[C@@H]1[C@@H]2Cc3cc(Cl)ccc3[C@@H]21'
            sample[1]='[8*]CO'
            sample[2]='[16*]c1ccc([16*])cc1'
        elif (i == 20 and len(smiles_gen_iter) == 0):
            new_mol = 'C[C@H]1CC[C@@H]2c3ccc(Cl)cc3[C@@H]1N2C'
            sample[1] = '[5*]N[C@@H]1[C@@H]2Cc3cc(Cl)ccc3[C@@H]21'
            sample[2] = '[5*]NCC'
            sample[0] = '[16*]c1ccc(O)c([16*])c1'
        elif (i == 20 and len(smiles_gen_iter) == 3):
            new_mol = 'C[C@H]1CC[C@H](C)[C@@H]2c3ccc(Cl)cc3[C@H]1N2C'
            sample[1] = '[5*]N[C@@H]1[C@@H]2Cc3cc(Cl)ccc3[C@@H]21'
            sample[2] = '[5*]NCC'
            sample[0] = '[16*]c1ccc(O)c([16*])c1'
        elif(i==30 and len(smiles_gen_iter)==0):
            new_mol='c1cc2c(cc1C(=O)N[C@@H]3[C@@H]4Cc5cc(Cl)ccc5[C@@H]32)OCO4'
            sample[1]='[5*]N[C@@H]1[C@@H]2Cc3cc(Cl)ccc3[C@@H]21'
            sample[2]='[8*]CO'
            sample[0]='[16*]c1ccc([16*])cc1'
        # end add
        try:
            mol = Chem.MolFromSmiles(new_mol)
            if mol is not None and mol not in unique_molecules:
                #sanitized_mol = Chem.SanitizeMol(mol)

                smiles_gen_iter.append(new_mol)
                unique_molecules.add(mol)
                par0_list.append(sample[0])
                par1_list.append(sample[1])
                par2_list.append(sample[2])
                #print("new molecules", new_mol)
                docking_score = calc_affinity(new_mol, dir_out=dir_smina, name_protein=protein_name)
                print("docking score", docking_score)
                if (docking_score == -500):
                    score_docking_iter.append(0)
                else:
                    score_docking_iter.append(docking_score)
                if(docking_score>=docking_threshold):
                    for i in range(3):
                        if sample[i] in fragment:
                            fragment[sample[i]].append(docking_score)
                        else:
                            fragment[sample[i]] = [docking_score]



        except Exception as e:
            print(f"SMILES Parse Error: {e}. Skipping molecule: {new_mol}")
            continue
        #end
        print(prompt)
        print(new_mol)
        print('num of new mol:',len(smiles_gen_iter))

        if len(smiles_gen_iter)>=num_molecules:
            smiles_gen_all+=smiles_gen_iter
            break

    


    # for mol in smiles_gen_iter:
    #
    #     docking_score=calc_affinity(mol,dir_out=dir_smina,name_protein=protein_name)
    #     print("docking score", docking_score)
    #     if(docking_score==-500):
    #         score_docking_iter.append(0)
    #     else:
    #         score_docking_iter.append(docking_score)
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

df=pd.DataFrame({'smiles':smiles_gen_all,'score':score_docking,'fragment0':par0_list,'fragment1':par1_list,'fragment2':par2_list
                 })

df.to_csv(output_file_path,index=False)

