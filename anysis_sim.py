from rdkit import DataStructs,Chem
from rdkit.Chem import AllChem
import numpy as np
from rdkit.Chem import MACCSkeys


from tools.tools import get_fingerpoint
# method=['MACCS','Topological','Morgan','Avalon']
method='MACCS'


def sim_tonimoto(fp_1,fp_2):
    common = 0
    num1=0
    num2=0
    # 判断有没有相同的数据, 没有相同数据则返回0
    for i in range(len(fp_1)):
        if fp_1[i]==1:
            num1=num1+1
        if fp_2[i]==1:
            num2=num2+1
        if ((fp_1[i]==1) and (fp_2[i]==1)):
            common=common+1



    common_num = common
    user1_num = num1
    user2_num = num2

    res = float(common_num) / (user1_num + user2_num - common_num)

    return res




def get_fp_list(smiles_list):
    fp_list = []
    for s in smiles_list:
        try:
            fp_list.append(get_fingerpoint())
        except Exception as e:
            print(f'Cont gen fp : {s}')
            continue
    return fp_list
def get_smiles_list(file_path):
    import pandas as pd



    a = pd.read_csv(file_path)
    smiles=list(a.iloc[:,0])
    return smiles
def fp_to_victor(fp_list):
    fp_victor = []
    for fp_tem in fp_list :
        fp_filted=[]
        for i in fp_tem:
            fp_filted.append(int(i))
        fp_victor.append(fp_filted)
    return fp_victor

if __name__ =='__main__':
    #每轮之间
    # file_list=['/home/th2024/lmlf/LMLF-main/4o_random_8UOB_100*20_3/init.csv','/home/th2024/lmlf/LMLF-main/4o_random_8UOB_100*20_3/result.csv']
    # fp=[]
    # label=[]
    # for i in range(len(file_list)):
    #     file_path_init = file_list[i]
    #     smiles_init = get_smiles_list(file_path_init)
    #     fp_init = get_fingerpoint(smiles_init, method)
    #     label_init = [i for j in fp_init]
    #     fp=fp+fp_init
    #     label=label+label_init
    # fp_victor=[]
    # for fp_tem in fp :
    #     fp_filted=[]
    #     for i in fp_tem:
    #         fp_filted.append(int(i))
    #     fp_victor.append(fp_filted)
    # fp=fp_victor
    # rear_begin=0
    # rear_end=100
    # last_begin=100
    # last_end=200
    # end=0
    # sim=0
    # num=0
    # while True:
    #     fp_rear=fp[rear_begin:rear_end]
    #     fp_last=fp[last_begin:last_end]
    #     for i in fp_rear:
    #         for j in fp_last:
    #             sim_tem=sim_tonimoto(i,j)
    #             sim=sim+sim_tem
    #             num=num+1
    #     rear_begin=rear_begin+100
    #     rear_end=rear_end+100
    #     last_begin=last_begin+100
    #     last_end=last_end+100
    #     if end==1:
    #         break
    #     if last_end>len(fp):
    #         last_end=len(fp)
    #     end=1
    # result=sim/num
    # print(result)
    #和初始
    file_init='/home/th2024/lmlf/LMLF-main/TSNE_8UOB_ours/top_500.csv'
    file_gen='/home/th2024/lmlf/LMLF-main/TSNE_8UOB_ours/top_500.csv'
    smiles_init = get_smiles_list(file_init)
    #smiles_init=smiles_init[:2000]
    fp_init = get_fingerpoint(smiles_init, method)
    fp_init=fp_to_victor(fp_init)
    smiles_gen = get_smiles_list(file_gen)
    #smiles_init = smiles_gen[:2000]
    fp_gen = get_fingerpoint(smiles_gen, method)
    fp_gen=fp_to_victor(fp_gen)
    num=0
    sim=0
    for fp1 in fp_init:
        for fp2 in fp_gen:
            sim_tem=sim_tonimoto(fp1,fp2)
            sim=sim+sim_tem
            num=num+1
    mean_sim=sim/num
    print(mean_sim)

