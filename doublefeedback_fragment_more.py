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
exp_path='fragment_PRMT5'

#old docking needed
directory_path = './LMLF/'
protein_name='7L1G'
num_molecules = 4
temperature = 0.3 # Controls the "creativity" of the generated molecules
num_generations = 100# The number of times to generate new molecules and feed them back into the modelcond
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
docking_threshold =6.5
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
        prompt = f'Generate a novel valid molecule SMILES which contains one fragment of [ {sample[0]} , {sample[1]},{sample[2]} ] at least and do not generate any English text,It is important to emphasize that a valid molecular formula does not contain free ions or any punctuations'

        import os
        from openai import OpenAI
        from openai.types.chat import completion_create_params

        os.environ["OPENAI_BASE_URL"] = "https://key.wenwen-ai.com/v1"
        os.environ["OPENAI_API_KEY"] = "sk-zY6KzDoXeI84S2xe95D0D4249b3e4a3b8d61F5800446F5Cd"
        model_engine = 'gpt-3.5-turbo'  # You can choose a different model if desired
        client = OpenAI()
        completion = client.chat.completions.create(
            model=model_engine,
            messages=[
                {"role": "user", "content": prompt}],
            n=1,
            max_tokens=60,
            temperature=0.6,
            stop="!",
            user="user"
        )

        new_mol=json.loads(completion.model_dump_json())["choices"][0]["message"]["content"]
        new_list = new_mol.split('.')
        for new in new_list:

            try:


                mol = Chem.MolFromSmiles(new)
                if mol is not None and mol not in unique_molecules:
                    #sanitized_mol = Chem.SanitizeMol(mol)

                    smiles_gen_iter.append(new)

                    par0_list.append(sample[0])
                    par1_list.append(sample[1])
                    par2_list.append(sample[2])
                    #print("new molecules", new_mol)
                    print(prompt)
                    print(new)
                    print('num of new mol:', len(smiles_gen_iter))

            except Exception as e:
                print(f"SMILES Parse Error: {e}. Skipping molecule: {new_mol}")
                continue
        #end


        if len(smiles_gen_iter)>=num_molecules:
            smiles_gen_all+=smiles_gen_iter
            break

    


    for mol in smiles_gen_iter:

        docking_score=calc_affinity(mol,dir_out=dir_smina,name_protein=protein_name)
        print("docking score", docking_score)
        if(docking_score==-500):
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

df=pd.DataFrame({'smiles':smiles_gen_all,'score':score_docking,'fragment0':par0_list,'fragment1':par1_list,'fragment2':par2_list
                 })

df.to_csv(output_file_path,index=False)

