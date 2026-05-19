import os
import subprocess
from loguru import logger
from rdkit import Chem
from rdkit.Chem import AllChem
import time
from easydict import EasyDict
from tqdm import tqdm
import numpy as np
#
#
# def calc_affinity(sml, name_protein='6GCT', dir_out='./', prefix='', os_type='linux', num_smina=1, cfg=None):
#
#     if cfg is not None:
#         autobox_add = cfg.autobox_add
#         seed = cfg.seed
#         exhaustiveness = cfg.exhaustiveness
#     else:
#         autobox_add = '16'
#         seed = '1000'
#         exhaustiveness = '16'
#
#     if name_protein == '6GCT':
#         file_protein = './datasets/{}_chainA_protein.pdbqt'.format(name_protein)
#         file_lig_ref = './datasets/{}_chainA_ligand.pdbqt'.format(name_protein)
#     else:
#         file_protein = './datasets/{}_chainA_protein.pdbqt'.format(name_protein)
#         file_lig_ref = './datasets/{}_chainA_ligand.pdbqt'.format(name_protein)
#
#     aff_array = np.zeros(num_smina).astype(np.float32)
#     for ii_num in range(num_smina):
#         try:
#             # smiles to pdb    'BOSearch/exp_BO_Search_seed-1/logs_smina/{}.pdb'.format（time）
#             mol = Chem.MolFromSmiles(sml)
#             m2 = mol
#             AllChem.EmbedMolecule(m2)
#             m3 = Chem.RemoveHs(m2)
#             file_output = os.path.join(dir_out, prefix + str(time.time()) + '.pdb')
#             Chem.MolToPDBFile(m3, file_output)
#
#             smina_cmd_output = os.path.join(dir_out, prefix + str(time.time()))
#             launch_args = ['/data/th2022/anaconda3/envs/lmlf/bin/smina', '-r', file_protein, '-l', file_output, '--autobox_ligand', file_lig_ref,
#                            '--autobox_add', autobox_add, '--seed', seed, '--exhaustiveness', exhaustiveness,
#                            '>>', smina_cmd_output]
#             #launch_args = ['smina', '-r', file_protein, '-l', file_output,
#             #                '--autobox_ligand', file_lig_ref, '--autobox_add', '10',
#             #                '--seed', '1000', '--exhaustiveness', '9', '-o', prefix+'dockres.pdb']
#             launch_string = ' '.join(launch_args)
#             logger.info(launch_string)
#             p = subprocess.Popen(launch_string, shell=True, stdout=subprocess.PIPE)
#
#             p.communicate()
#
#
#
#
#             affinity = 500
#             with open(smina_cmd_output, 'r') as f:
#                 for lines in f.readlines():
#                     lines = lines.split()
#                     if len(lines) == 4 and lines[0] == '1':
#                         affinity = float(lines[1])
#
#             if 'win' in os_type:
#                 prefix_del = 'del '
#             else:
#                 prefix_del = 'rm -rf '
#             p = subprocess.Popen(prefix_del + smina_cmd_output, shell=True, stdout=subprocess.PIPE)
#             p.communicate()
#             p = subprocess.Popen(prefix_del + file_output, shell=True, stdout=subprocess.PIPE)
#             p.communicate()
#         except:
#             affinity = 500
#
#         if affinity == 500:
#             logger.error('**** Affinity error ... ****')
#             aff_array = np.array([affinity] * num_smina).astype(np.float32)
#             break
#         aff_array[ii_num] = affinity
#
#     return - aff_array.mean()
#
#
#
def calc_affinity(mol_file,
                  name_protein='8UOB',
                  dir_out='./', prefix='', os_type='linux',
                  dock_type=0, cfg=None):

    file_protein='datasets/{}_chainA_protein.pdb'.format(name_protein)
    file_lig_ref='datasets/{}_chainA_ligand.pdb'.format(name_protein)
    if not os.path.exists(file_lig_ref):
        file_lig_ref = 'datasets/{}_chainA_ligand.sdf'.format(name_protein)

    flexdist='6'
    if cfg is not None:
        # flexdist = str(cfg['flexdist'])
        autobox_add = str(cfg['autobox_add'])
        seed = str(cfg['seed'])
        exhaustiveness = str(cfg['exhaustiveness'])
    else:
        flexdist = '6'
        autobox_add = '1'
        seed = '1000'
        exhaustiveness = '16'

    try:
        mol = Chem.MolFromPDBFile(mol_file)
        m2 = Chem.AddHs(mol)
        status=AllChem.EmbedMolecule(m2)
        if status == -1:
            raise ValueError('RDKIT fail')
        m3 = Chem.RemoveHs(m2)
        file_output = os.path.join(dir_out, prefix + str(time.time()) + '.pdb')
        Chem.MolToPDBFile(m3, file_output)


        smina_cmd_output = os.path.join(dir_out, prefix + str(time.time()))

        if dock_type == 0:
            launch_args = ['/home/th2024/anaconda3/envs/lmlf/bin/smina', '-r', file_protein, '-l', file_output, '--autobox_ligand', file_lig_ref,
                           '--autobox_add', autobox_add, '--seed', seed, '--exhaustiveness', exhaustiveness,
                           '-o', os.path.join(dir_out,prefix + 'smina_out.mol2'), '>>', smina_cmd_output]
        else:
            launch_args = ['/home/th2024/anaconda3/envs/lmlf/bin/smina', '-r', file_protein, '--flexdist_ligand', file_lig_ref, '--flexdist', flexdist,
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

        if 'win' in os_type:
            prefix_del = 'del '
        else:
            prefix_del = 'rm -rf '
        p = subprocess.Popen(prefix_del + smina_cmd_output, shell=True, stdout=subprocess.PIPE)
        p.communicate()
        p = subprocess.Popen(prefix_del + file_output, shell=True, stdout=subprocess.PIPE)
        p.communicate()
    except Exception as e:
        affinity = 500
        print(e)

    if affinity == 500:
        logger.error('**** Affinity error. ****')

    return -affinity

if __name__=='__main__':
    score=calc_affinity(mol_file='/home/th2024/lmlf/LMLF-main/datasets/7MJS_chainA_ligand.pdb',name_protein='7MJS')
    print(score)


