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


model_engine = 'xxx' # You can choose a different model if desired
temperature = 0.3 # Controls the "creativity" of the generated molecules
num_generations = 100# The number of times to generate new molecules and feed them back into the modelcond
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
from tools.tools import get_fragment
exp_path='PRMT5_fragment'

#old docking needed
directory_path = './LMLF/'
protein_name='7L1G'
# target_file_path = 'one-box/drd2.jsonl'
target_file_path =os.path.join(exp_path,'init.jsonl')
result_file_path=os.path.join(exp_path,'result.jsonl')
output_file_path = os.path.join(exp_path,'result.csv')
parent_file_path = os.path.join(exp_path,'parent.csv')
from tools.tools import makedir
makedir("smina1")
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
gen_list=[]
par0_list=[]
par1_list=[]
par2_list=[]
threshold_increment_frequency = 10
data = []
docking_threshold =8
unique_molecules = set()
fragment={}
# Generate and feed back new molecules k times
for i in range(1, num_generations):
    if i % threshold_increment_frequency == 0 and i > 0:
        docking_threshold += 1
    
    #Generating molecules
    new_molecules = []
    num_molecules = 4
    print("iteration", i)
    smiles_gen_iter=[]
    if (i==1):
        with jsonlines.open(target_file_path) as reader:
            for line in reader:
                if "\n" not in line:
                    if(line['label']=='1'):
                        target_molecules.append(line['smiles'])
                        target_labels.append(line['label'])
        new_molecules = target_molecules
    else:
        for point in new_target_molecules:
            target_molecules.append(point['smiles'])
            target_labels.append(point['label'])
        new_molecules=target_molecules
    print('num example:',len(target_labels))
    for molecule in new_molecules:
        try:
            score=calc_affinity(molecule,dir_out="smina",name_protein=protein_name)
            fragment_tem=get_fragment(molecule,'BRICS')
            for f in fragment_tem:
                if f in fragment :
                    fragment[f].append(score)
                else:
                    fragment[f]=[score]
        except Exception as e:
            print('fail to generate fragment')

    fragment_avg={}
    for key,valus in fragment.items():
        fragment_avg[key]=statistics.mean(valus)
    # fragment_avg_sorted=sorted(fragment_avg.items(),key=lambda d: d[1],reverse=True)
    population=np.array(list(fragment_avg.keys()))
    weights=np.array(list(fragment_avg.values()))

    while True:
        sample = np.random.choice(population, size=3, replace=False, p=weights / sum(weights))
        #target_mol = None

        # with jsonlines.open(target_file_path) as reader:
        #     for line in reader:
        #         target_mol = line['smiles']
        #         target_label = line['label']
        #         break
        print("avg dict:",fragment_avg)
        print("fragment choised:", sample)
        #change
        prompt = f'Generate a novel valid molecule SMILES which contains one fragments of [ {sample[0]} , {sample[1]},{sample[2]} ]at least and do not generate any English text'

        import os
        from openai import OpenAI
        from openai.types.chat import completion_create_params

        os.environ["OPENAI_BASE_URL"] = "https://key.wenwen-ai.com/v1"
        os.environ["OPENAI_API_KEY"] = "sk-R2iz32ySqmWIolon294bE059Ec974d3bAbC8E8F6760d8dB7"
        model_engine = 'gpt-3.5-turbo-instruct'  # You can choose a different model if desired
        client = OpenAI()
        completion = client.chat.completions.create(
            model='gpt-3.5-turbo',
            messages=[
                {"role": "user", "content": prompt}],
            n=1,
            max_tokens=60,
            temperature=0.6,
            stop="!",
            user="user"
        )

        new_mol=json.loads(completion.model_dump_json())["choices"][0]["message"]["content"]
        try:
            mol = Chem.MolFromSmiles(new_mol)
            if mol is not None and mol not in unique_molecules:
                #sanitized_mol = Chem.SanitizeMol(mol)
                smiles_gen_iter.append(new_mol)
                gen_list.append(new_mol)
                par0_list.append(sample[0])
                par1_list.append(sample[1])
                par2_list.append(sample[2])
                #print("new molecules", new_mol)
                print("new_molecules", gen_list)
        except Exception as e:
            print(f"SMILES Parse Error: {e}. Skipping molecule: {new_mol}")
            continue
        #end
        print('num of new mol:',len(smiles_gen_iter))
        if len(smiles_gen_iter)==num_molecules:
            new_molecules+=smiles_gen_iter
            break

    
    #clculate docking scores
    docking_scores = []

    for mol in new_molecules:
        #old
        #docking_score = calculate_docking_score(mol)
        docking_score=calc_affinity(mol,dir_out="smina",name_protein=protein_name)
        print("docking score", docking_score)
        if(docking_score==-500):
            docking_scores.append(None)
        else:
            docking_scores.append(docking_score)
    #calculate QED scores
    qed_scores = []
    for mol in new_molecules:
        qed_score = calculate_qed_score(mol)
        print(mol,"QED score", qed_score)
        qed_scores.append(qed_score)
    
    #calculate RADcores
    mw_scores = []
    logp_scores = []
    radscores = []
    for mol in new_molecules:
        # try:
        mw = Descriptors.MolWt(Chem.MolFromSmiles(mol))
        logp = Descriptors.MolLogP(Chem.MolFromSmiles(mol))
        # sas = rdmd.SyntheticAccessibility(mol)
        mw_scores.append(mw)
        logp_scores.append(logp)
        radscore = calculate_radscore(mol)
        print("RD score", radscore)
        radscores.append(radscore)
        # except:
        #     continue
    
    #select molecules based on upper quartiles of docking scores and RADcores
    filtered_docking_scores = [score for score in docking_scores if score is not None]
    filtered_mw_scores = [score for score in mw_scores if score is not None]
    filtered_logp_scores = [score for score in logp_scores if score is not None]
    filtered_qed_scores = [score for score in qed_scores if score is not None]
    print("radscores", radscores)
    filtered_rad_scores = [score for score in radscores if score is not None] 
    #docking_threshold = np.percentile(filtered_docking_scores, 75)
    
    print("docking threshold", docking_threshold)
    #qed_threshold = np.percentile(filtered_qed_scores, 75)
    qed_threshold = 4
    if len(filtered_rad_scores) > 0:
        radscore_threshold = np.percentile(filtered_rad_scores, 75)
    else:
        radscore_threshold = 0.0  # Assign a default value or handle the case appropriately
    #radscore_threshold = np.percentile(filtered_rad_scores, 75)
    labels = []
    mw_threshold = 700
    logp_threshold = 6.0
    new_target_molecules = []
    for mol, docking_score, qed_score, radscore,mw ,logp in zip(new_molecules, filtered_docking_scores, filtered_qed_scores, filtered_rad_scores,mw_scores,logp_scores):
        if (
            docking_score >= docking_threshold 
            and radscore >= radscore_threshold
            and mw <= mw_threshold
            and logp <= logp_threshold
        ):
            labels.append('1')
            new_target_molecules.append({'smiles': mol, 'label': '1'})
        else:
            labels.append('0')
    
    #Repeat or stop based on joint median change
    # if i >= 1:
    #     prev_joint_median = np.median([max(docking_score, qed_score, radscore) for docking_score, qed_score, radscore in zip(filtered_docking_scores, filtered_qed_scores, filtered_rad_scores)])
    #
    #     joint_median = np.median(filtered_docking_scores + filtered_rad_scores + filtered_qed_scores)
    #     if abs(joint_median - prev_joint_median) < 0.01:
    #         break
    

    data.extend(list(zip(new_molecules, filtered_docking_scores, filtered_qed_scores, filtered_rad_scores, labels)))
    # Perform further operations or analysis with the data as needed
    print("data", data)


    #with jsonlines.open(target_file_path, mode='a') as writer:
        # writer.write('\\n')
        # writer.write_all(new_target_molecules)
    with jsonlines.open(result_file_path, mode='a') as writer:
        # writer.write("\\n")
        for molecule in new_target_molecules:
            writer.write(molecule)
            writer.write('\n')
    # with open(target_file_path, mode='a') as outfile:
    #     for hostDict in target_file_path:
    #         json.dump(hostDict, outfile)
    #         outfile.write('\n')

df = pd.DataFrame(data, columns=['Molecule', 'Docking Score', 'QED Score', 'RADscore', 'Label'])
df.to_csv(output_file_path, index=False)
joint_scores = filtered_docking_scores + filtered_rad_scores + filtered_qed_scores
pd.DataFrame({'gen_chain':gen_chian}).to_csv(os.path.join(exp_path, 'gen_chain.csv'),index=False)
# Print and analyze the results
# print(f'generation {i}:')
# print(f'joint median: {joint_median}')
# print('generated molecules:')
# print(data)
pd.DataFrame({'gen_smiles':gen_list,'fragment0':par0_list,'fragment1':par1_list,'fragment2':par2_list}).to_csv(parent_file_path,index=False)