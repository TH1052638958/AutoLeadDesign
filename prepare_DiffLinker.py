import os.path


from tools.tools import get_fragment
import pandas as pd
import os
import subprocess
from loguru import logger
from rdkit import Chem
from rdkit.Chem import AllChem
import time
from tools.tools import makedir
def calc_affinity(sml,
                  name_protein='8UOB',
                  dir_out='./', prefix='', os_type='linux',
                  dock_type=0, cfg=None,mol2file='1.mol2'):

    file_protein='datasets/{}_chainA_protein.pdb'.format(name_protein)
    file_lig_ref='datasets/{}_chainA_ligand.pdb'.format(name_protein)


    if cfg is not None:
        # flexdist = str(cfg['flexdist'])
        autobox_add = str(cfg['autobox_add'])
        seed = str(cfg['seed'])
        exhaustiveness = str(cfg['exhaustiveness'])
    else:
        flexdist = '6'
        autobox_add = '6'
        seed = '1000'
        exhaustiveness = '16'

    try:
        mol = Chem.MolFromSmiles(sml)
        m2 = Chem.AddHs(mol)
        AllChem.EmbedMolecule(m2)
        m3 = Chem.RemoveHs(m2)
        file_output = os.path.join(dir_out, prefix + str(time.time()) + '.pdb')
        Chem.MolToPDBFile(m3, file_output)

        smina_cmd_output = os.path.join(dir_out, prefix + str(time.time()))


        if dock_type == 0:
            launch_args = ['/home/th2024/anaconda3/envs/lmlf/bin/smina', '-r', file_protein, '-l', file_output, '--autobox_ligand', file_lig_ref,
                           '--autobox_add', autobox_add, '--seed', seed, '--exhaustiveness', exhaustiveness,
                           '-o', mol2file, '>>', smina_cmd_output]
        else:
            launch_args = ['smina', '-r', file_protein, '--flexdist_ligand', file_lig_ref, '--flexdist', flexdist,
                           '-l', file_output, '--autobox_ligand', file_lig_ref, '--autobox_add', autobox_add,
                           '--seed', seed, '--exhaustiveness', exhaustiveness,
                           '-o', os.path.join(dir_out,prefix + 'smina_out.mol2'), '>>', smina_cmd_output]

        # launch_args = ['smina', '-r', file_protein, '-l', file_output,
        #                '--autobox_ligand', file_lig_ref, '--autobox_add', '10',
        #                '--seed', '1000', '--exhaustiveness', '9', '-o', prefix+'dockres.pdb']
        launch_string = ' '.join(launch_args)
        logger.info(launch_string)
        p = subprocess.Popen(launch_string, shell=True, stdout=subprocess.PIPE)
        p.communicate()

        affinity = 500
        with open(smina_cmd_output, 'r') as f:
            for lines in f.readlines():
                lines = lines.split()
                if len(lines) == 4 and lines[0] == '1':
                    affinity = float(lines[1])
        molecules = Chem.MolFromMolFile(mol2file)
        first_conformer = molecules.GetConformer(0)
        first_molecule = Chem.Mol(molecules)
        first_molecule.RemoveAllConformers()
        first_molecule.AddConformer(first_conformer)
        if 'win' in os_type:
            prefix_del = 'del '
        else:
            prefix_del = 'rm -rf '
        p = subprocess.Popen(prefix_del + smina_cmd_output, shell=True, stdout=subprocess.PIPE)
        p.communicate()
        p = subprocess.Popen(prefix_del + mol2file, shell=True, stdout=subprocess.PIPE)
        p.communicate()
        p = subprocess.Popen(prefix_del + file_output, shell=True, stdout=subprocess.PIPE)
        p.communicate()
    except:
        affinity = 500

    if affinity == 500:
        logger.error('**** Affinity error. ****')
        affinity=0
    else:
        output_file = mol2file
        Chem.MolToMolFile(first_molecule, output_file)
    return -affinity




import argparse

parser = argparse.ArgumentParser(description='preDiffLinker')
parser.add_argument('--exp-path', type=str, default='DL_8UOB')
parser.add_argument('--cfg-smina', type=str, default='config_smina_8UOB.yaml')
parser.add_argument('--protein-name', type=str, default='8UOB')
parser.add_argument('--fragment-process', type=bool, default=True)
parser.add_argument('--fragment', type=str, default='BRICS')
args = parser.parse_args()
import ruamel.yaml as yaml
from easydict import EasyDict
dir_smina=os.path.join(args.exp_path,'smina')
makedir(dir_smina)
with open(args.cfg_smina) as f:
    cfg_smina = yaml.load(f, Loader=yaml.Loader)
args.cfg_smina = EasyDict(cfg_smina)
protein_name=args.protein_name
Smiles_path=os.path.join(args.exp_path,'init.csv')
data=df_Smiles=pd.read_csv(Smiles_path)
smiles=list(data.iloc[:,0])
fragment=[]
score=[]
for s in smiles:
    fragment_tem = get_fragment(s, args.fragment, args.fragment_process)
    for f in fragment_tem:
        if (f in fragment):
            continue
        else:
            fragment.append(f)
for n in range(len(fragment)):
    f=fragment[n]
    outfile=os.path.join(args.exp_path,'fragment_{}.mol'.format(str(n)))
    aff=calc_affinity(f,protein_name,dir_smina,mol2file=outfile)
    score.append(aff)

result_file=os.path.join(args.exp_path,'fragment_score.csv')
pd.DataFrame({'fragment':fragment,'score':score}).to_csv(result_file,index=False)

