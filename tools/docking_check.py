import os
import subprocess
from loguru import logger
from rdkit import Chem
from rdkit.Chem import AllChem
import time


def calc_affinity(sml,
                  file_protein='./test_pdbs/1a9u/1a9u_protein.pdb',
                  file_lig_ref='./test_pdbs/1a9u/1a9u_ligand.sdf',
                  dir_out='./', prefix='', os_type='linux',
                  dock_type=0, cfg=None):

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
                           '-o', prefix + 'smina_out.mol2', '>>', smina_cmd_output]
        else:
            launch_args = ['smina', '-r', file_protein, '--flexdist_ligand', file_lig_ref, '--flexdist', flexdist,
                           '-l', file_output, '--autobox_ligand', file_lig_ref, '--autobox_add', autobox_add,
                           '--seed', seed, '--exhaustiveness', exhaustiveness,
                           '-o', prefix + 'smina_out.mol2', '>>', smina_cmd_output]

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
