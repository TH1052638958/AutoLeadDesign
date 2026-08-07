import os
import pandas as pd
import argparse
import logging


def makedir(path):
    try:
        os.makedirs(path)
    except OSError:
        pass


def load_dataset(dir_data='./datasets', name_database='ZINC'):
    if 'ZINC' in name_database:
        dataset = pd.read_csv(os.path.join(dir_data, '{}.txt'.format(name_database)), sep='\t')
        # return dataset.values
        return dataset
    else:
        print('Now, the database of {} is not supported ...'.format(name_database))
        exit(1)


def init_args_data(dataset):
    """
    Reads in the arguments for the script for a given dataset.

    Parameters
    ----------
    dataset : :class:`str`
        Dataset being used.  Currently 'ZINC' is supported.

    Returns
    -------
    args : :class:`Namespace`
        Namespace with a dictionary of arguments where the key is the name of
        the argument and the item is the input value.
    """

    parser = argparse.ArgumentParser(description='Dataset'.format(dataset))
    parser.add_argument('--name-data', type=str, default='ZINC', help='dataset name')
    args_data = parser.parse_args()

    return args_data


class Logger(object):
    def __init__(self, log_file_name, log_level, logger_name):
        self.__logger = logging.getLogger(logger_name)
        self.__logger.setLevel(log_level)
        file_handler = logging.FileHandler(log_file_name)
        console_handler = logging.StreamHandler()
        # formatter = logging.Formatter(
        #     '[%(asctime)s] - [%(filename)s line:%(lineno)d] : %(message)s')
        formatter = logging.Formatter('[%(asctime)s] - : %(message)s')
        file_handler.setFormatter(formatter)
        console_handler.setFormatter(formatter)
        self.__logger.addHandler(file_handler)
        self.__logger.addHandler(console_handler)

    def get_log(self):
        return self.__logger
