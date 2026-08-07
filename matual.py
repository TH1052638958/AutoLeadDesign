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


model_engine = 'gpt-3.5-turbo' # You can choose a different model if desired
temperature = 0.3 # Controls the "creativity" of the generated molecules
num_generations = 20# The number of times to generate new molecules and feed them back into the modelcond
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
from rdkit.Contrib.SA_Score import sascorer
#change to pre line
# sys.path.append(os.path.join(RDConfig.RDContribDir, 'SA_Score'))
# import sascorer
#change end
from tools.docking import calc_affinity
exp_path='PRMT5_1'

#old docking needed
directory_path = './LMLF/'
# target_file_path = 'one-box/drd2.jsonl'
target_file_path =os.path.join(exp_path,'init.jsonl')
result_file_path=os.path.join(exp_path,'result.jsonl')
output_file_path = os.path.join(exp_path,'result.csv')
from tools.tools import makedir
makedir("smina")
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
            

      



# target_mol = None
# with jsonlines.open(target_file_path) as reader:
#     for line in reader:
#         target_mol = line['smiles']
#         target_label = line['label']
        # break
target_molecules = []
target_labels = []
gen_chian=[]
threshold_increment_frequency = 10
data = []
docking_threshold =9
unique_molecules = set()
target_mol='C[C@H](CNC(=O)c1cccc(-c2ccoc2)c1)CC(=O)O'
target_label=1
prompt = f'Generate a novel valid molecule similar to {target_mol} that is {target_label}-class and do not generate any English text'

import os
from openai import OpenAI


os.environ["OPENAI_BASE_URL"] = "https://key.wenwen-ai.com/v1"
os.environ["OPENAI_API_KEY"] = "sk-EF2D4jHDHF9psLq023395862Fd8c4a8d83Dd428b597dEbDc"
model_engine = 'gpt-3.5-turbo-instruct'  # You can choose a different model if desired
client = OpenAI()
completion = client.chat.completions.create(
            model='gpt-3.5-turbo',
            messages=[
                {"role": "system", "content": prompt}],
            n=1,
            max_tokens=60,
            temperature=0.6,
            stop="!",
            user="test"
        )

new_mol=json.loads(completion.model_dump_json())["choices"][0]["message"]["content"]
print('new mol:',new_mol)
docking_score=calc_affinity(new_mol,dir_out="smina")
print("docking score", docking_score)
