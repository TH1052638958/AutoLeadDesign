import os
import sys
import time
import argparse
import logging
import random
import pandas as pd
import numpy as np
from tqdm import tqdm

import utils
from tools import docking_check as docking


parser = argparse.ArgumentParser(description='ZINC Docking')
parser.add_argument('--dir-exp', type=str, default='exp_Check', help='dir for experiment')
parser.add_argument('--dir-data', type=str, default='data_8uob', help='dir for dataset')
parser.add_argument('--name-protein', type=str, default='8uob', help='name of protein (default: 1a9u, 6GCT, 1iep)')
parser.add_argument('--name-sml', type=str, default='8uob', help='name of the database to search')
parser.add_argument('--seed', type=int, default=1, help='name of the database to search')
parser.add_argument('--dock-type', type=int, default=0,
                    help='the type of docking (default: 0), 0 for not flex, 1 for flex')

args = parser.parse_args()

# use parser to get the training information
if __name__ == "__main__":
    dir_root = os.getcwd()

    args.name_sml = args.name_sml + '-{}'.format(args.seed)
    args.os = sys.platform
    args.dir_exp = os.path.join(dir_root, args.dir_exp + '-' + args.name_protein)
    args.dir_data = os.path.join(dir_root, args.dir_data)
    args.dir_smina = os.path.join(args.dir_exp, 'logs_smina')
    # creating work dirs and initializing some variables
    utils.makedir(args.dir_exp)
    utils.makedir(args.dir_smina)

    seeds = [1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000, 10000]
    cfg = {'flexdist': 6, 'autobox_add': 16, 'exhaustiveness': 16, 'seed': 1000}

    timestamp = time.strftime('%Y%m%d_%H%M%S', time.localtime())
    log_file = os.path.join(args.dir_exp, '{}.log'.format(timestamp))
    logger = utils.Logger(log_file_name=log_file, log_level=logging.DEBUG, logger_name='SMINA').get_log()
    logger.info('---- Args Info ----')
    logger.info(args)
    logger.info('---- Smina Config ----')
    logger.info(cfg)

    # load dataset
    data = utils.init_args_data(args.name_sml)
    if '6GCT' in args.name_protein:
        data.protein = os.path.join(args.dir_data, '{}_chainA_protein.pdbqt'.format(args.name_protein))
        data.ligand = os.path.join(args.dir_data, '{}_chainA_ligand.pdbqt'.format(args.name_protein))
    elif '1iep' in args.name_protein:
        data.protein = os.path.join(args.dir_data, '{}_chainA.pdbqt'.format(args.name_protein))
        data.ligand = os.path.join(args.dir_data, '{}_chainA_ligand.pdbqt'.format(args.name_protein))
    else:
        data.protein = os.path.join(args.dir_data, '{}_protein.pdbqt'.format(args.name_protein))
        data.ligand = os.path.join(args.dir_data, '{}_ligand.pdbqt'.format(args.name_protein))
    data.sml_data_all = pd.read_csv(os.path.join(args.dir_data, '{}.csv'.format(args.name_sml)))
    data.sml_data = data.sml_data_all.values[:, 1]
    data.num_data = len(data.sml_data)
    data.index_data = list(range(0, data.num_data))

    num_time = 1
    aff_array = np.zeros([data.num_data, num_time])
    for ii_time in range(num_time):
        # cfg['seed'] = seeds[ii_time]
        start_time = time.time()
        pbar = tqdm(range(data.num_data), total=data.num_data)
        for ii_idx in pbar:
            pbar.set_description('processing -> time: {}/{}'.format(ii_time+1, num_time))
            aff_array[ii_idx, ii_time] = docking.calc_affinity(data.sml_data[ii_idx],
                                                               file_protein=data.protein,
                                                               file_lig_ref=data.ligand,
                                                               dir_out=args.dir_smina,
                                                               prefix='exp-add{}-{}_sml-{}_'.format(cfg['autobox_add'],
                                                                                                    args.seed, ii_idx),
                                                               os_type=args.os,
                                                               dock_type=args.dock_type,
                                                               cfg=cfg)
        end_time = time.time()
        cost_time = end_time - start_time
        logger.info('calc_affinity at {} time, total time cost is {:.6f}, '
                    'the mean time for each docking is {:.6f}'.format(ii_time+1, cost_time, cost_time/data.num_data))
        # data.sml_data.insert(loc=1, column='aff_{}'.format(ii_time+1), value=aff_array[:, ii_time])
        # data.sml_data.to_csv(os.path.join(args.dir_exp, '{}-{}.txt'.format(args.name_database, args.name_protein)),
        #                      sep='\t', index=False)
    aff_mean = np.mean(aff_array, axis=1)
    data.sml_data_all.insert(loc=2, column='aff_mean', value=aff_mean)
    data.sml_data_all.to_csv(os.path.join(args.dir_exp, '{}_{}.csv'.format(args.name_sml, args.name_protein)),
                             sep='\t', index=False)
