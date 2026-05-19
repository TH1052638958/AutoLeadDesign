import os
import subprocess

import rdkit.Chem.rdmolfiles
from loguru import logger
from rdkit import Chem
from rdkit.Chem import AllChem
import time
from tools import docking_check
from tools import docking


# def calc_affinity(sml,
#                   file_protein='./test_pdbs/1a9u/1a9u_protein.pdb',
#                   file_lig_ref='./test_pdbs/1a9u/1a9u_ligand.sdf',
#                   dir_out='./', prefix='', os_type='linux',
#                   dock_type=0, cfg=None):
#
#     if cfg is not None:
#         # flexdist =0
#         autobox_add = str(cfg['autobox_add'])
#         seed = str(cfg['seed'])
#         exhaustiveness = str(cfg['exhaustiveness'])
#     else:
#         flexdist = '6'
#         autobox_add = '6'
#         seed = '1000'
#         exhaustiveness = '16'
#
#     try:
#         mol = Chem.MolFromSmiles(sml)
#         m2 = Chem.AddHs(mol)
#         AllChem.EmbedMolecule(m2)
#         m3 = Chem.RemoveHs(m2)
#         file_output = os.path.join(dir_out, prefix + str(time.time()) + '.pdb')
#         Chem.MolToPDBFile(m3, file_output)
#
#         smina_cmd_output = os.path.join(dir_out, prefix + str(time.time()))
#
#         if dock_type == 0:
#             launch_args = ['/home/th2024/anaconda3/envs/lmlf/bin/smina', '-r', file_protein, '-l', file_output, '--autobox_ligand', file_lig_ref,
#                            '--autobox_add', autobox_add, '--seed', seed, '--exhaustiveness', exhaustiveness,
#                            '-o', prefix + 'smina_out.mol2', '>>', smina_cmd_output]
#         else:
#             launch_args = ['smina', '-r', file_protein, '--flexdist_ligand', file_lig_ref, '--flexdist', flexdist,
#                            '-l', file_output, '--autobox_ligand', file_lig_ref, '--autobox_add', autobox_add,
#                            '--seed', seed, '--exhaustiveness', exhaustiveness,
#                            '-o', prefix + 'smina_out.mol2', '>>', smina_cmd_output]
#
#         # launch_args = ['smina', '-r', file_protein, '-l', file_output,
#         #                '--autobox_ligand', file_lig_ref, '--autobox_add', '10',
#         #                '--seed', '1000', '--exhaustiveness', '9', '-o', prefix+'dockres.pdb']
#         launch_string = ' '.join(launch_args)
#         logger.info(launch_string)
#         p = subprocess.Popen(launch_string, shell=True, stdout=subprocess.PIPE)
#         p.communicate()
#
#         affinity = 500
#         with open(smina_cmd_output, 'r') as f:
#             for lines in f.readlines():
#                 lines = lines.split()
#                 if len(lines) == 4 and lines[0] == '1':
#                     affinity = float(lines[1])
#
#
#     except:
#         affinity = 500
#
#     if affinity == 500:
#         logger.error('**** Affinity error. ****')
#
#     return -affinity
# def calc_affinity(sml,
#                   idx='0',
#                   file_protein='./test_pdbs/1a9u/1a9u_protein.pdb',
#                   file_lig_ref='./test_pdbs/1a9u/1a9u_ligand.sdf',
#                   dir_out='./', prefix='', os_type='linux',
#                   dock_type=0, cfg=None):
#
#     if cfg is not None:
#         # flexdist = str(cfg['flexdist'])
#         autobox_add = str(cfg['autobox_add'])
#         seed = str(cfg['seed'])
#         exhaustiveness = str(cfg['exhaustiveness'])
#     else:
#         flexdist = '6'
#         autobox_add = '6'
#         seed = '1000'
#         exhaustiveness = '16'
#
#     try:
#         mol = Chem.MolFromSmiles(sml)
#
#         m2 = Chem.AddHs(mol)
#         AllChem.EmbedMolecule(m2,randomSeed=10)
#
#         m3 = Chem.RemoveHs(m2)
#
#         file_output = os.path.join(dir_out, prefix + str(time.time()) + '.pdb')
#         Chem.MolToPDBFile(m3, file_output)
#
#         smina_cmd_output = os.path.join(dir_out, prefix + str(time.time()))
#
#         if dock_type == 0:
#             launch_args = ['/home/th2024/anaconda3/envs/lmlf/bin/smina', '-r', file_protein, '-l', file_output, '--autobox_ligand', file_lig_ref,
#                            '--autobox_add', autobox_add, '--seed', seed, '--exhaustiveness', exhaustiveness,
#                            '-o', os.path.join(dir_out,prefix + f'smina_out_{idx}.pdbqt'), '>>', smina_cmd_output]
#         else:
#             launch_args = ['smina', '-r', file_protein, '--flexdist_ligand', file_lig_ref, '--flexdist', flexdist,
#                            '-l', file_output, '--autobox_ligand', file_lig_ref, '--autobox_add', autobox_add,
#                            '--seed', seed, '--exhaustiveness', exhaustiveness,
#                            '-o', prefix + 'smina_out.mol2', '>>', smina_cmd_output]
#
#         # launch_args = ['smina', '-r', file_protein, '-l', file_output,
#         #                '--autobox_ligand', file_lig_ref, '--autobox_add', '10',
#         #                '--seed', '1000', '--exhaustiveness', '9', '-o', prefix+'dockres.pdb']
#         launch_string = ' '.join(launch_args)
#         logger.info(launch_string)
#         p = subprocess.Popen(launch_string, shell=True, stdout=subprocess.PIPE)
#         p.communicate()
#
#         affinity = 500
#         with open(smina_cmd_output, 'r') as f:
#             for lines in f.readlines():
#                 lines = lines.split()
#                 if len(lines) == 4 and lines[0] == '1':
#                     affinity = float(lines[1])
#
#         # if 'win' in os_type:
#         #     prefix_del = 'del '
#         # else:
#         #     prefix_del = 'rm -rf '
#         # p = subprocess.Popen(prefix_del + smina_cmd_output, shell=True, stdout=subprocess.PIPE)
#         # p.communicate()
#         # p = subprocess.Popen(prefix_del + file_output, shell=True, stdout=subprocess.PIPE)
#         # p.communicate()
#     except:
#         affinity = 500
#
#     if affinity == 500:
#         logger.error('**** Affinity error. ****')
#
#     return -affinity
def calc_affinity(sml,
                  name_protein='8UOB',
                  dir_out='./', prefix='', os_type='linux',
                  dock_type=0, cfg=None):

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
                           '-o', os.path.join(dir_out,prefix + f'smina_out_{str(time.time())}.mol2'), '>>', smina_cmd_output]
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

        if 'win' in os_type:
            prefix_del = 'del '
        else:
            prefix_del = 'rm -rf '
        p = subprocess.Popen(prefix_del + smina_cmd_output, shell=True, stdout=subprocess.PIPE)
        p.communicate()
        p = subprocess.Popen(prefix_del + file_output, shell=True, stdout=subprocess.PIPE)
        p.communicate()
    except:
        affinity = 500

    if affinity == 500:
        logger.error('**** Affinity error. ****')

    return -affinity


if __name__=='__main__':
    #CC1=CC=CN=C1N[C@@H]1CCC2CCCCC21N=Cc1ccc2c(c1)CNC2c1nc(N)c2ccccc2n1
    smiles='C[C@H]1CC[C@H](C2=C3C=CC=CC=C3NC=C2)C[C@H]1[C@H]1CN(C2=C3C=CC=CC=C3NC=C2)C[C@@H]1C'
    mol = Chem.MolFromSmiles(smiles)


    import ruamel.yaml as yaml
    from easydict import EasyDict
    with open('config_smina_8UOB.yaml') as f:
        cfg_smina = yaml.load(f, Loader=yaml.Loader)
    cfg_smina = EasyDict(cfg_smina)
    cfg_smina['autobox_add']=1
    # a = docking_check.calc_affinity(smiles,
    #                                 file_protein='datasets/8UOB_chainA_protein.pdb',
    #                                 file_lig_ref='datasets/8UOB_chainA_ligand.pdb',
    #                                 dir_out='smina',
    #
    #                                 cfg=cfg_smina)
    # print(a)
    for i in range(5):
        # a=calc_affinity(smiles,
        #                       file_protein='datasets/8UOB_chainA_protein.pdb',
        #                       file_lig_ref='datasets/8UOB_chainA_ligand.pdb',
        #                       dir_out='smina',idx=str(i),
        #
        #
        #
        #
        #                       cfg=cfg_smina)
        # print(a)
        print(calc_affinity(smiles, dir_out='smina', name_protein='8UOB', cfg=cfg_smina))
        # print(docking.calc_affinity(smiles,dir_out='smina',name_protein='8UOB',cfg=cfg_smina))

        print('-------------------------------------------------------')
