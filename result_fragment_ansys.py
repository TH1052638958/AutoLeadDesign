import pandas as pd

from tools.tools import makedir
import argparse
import os
import numpy as np
from rdkit.Chem import Draw
from PIL import Image
from rdkit import Chem







parser = argparse.ArgumentParser(description='LMLF')
parser.add_argument('--result-dir', type=str, default='random_8UOB_fragment')
parser.add_argument('--iter', type=int, default=40)
parser.add_argument('--num-fragment', type=int, default=20)
parser.add_argument('--output-dir', type=str, default='fragment_8UOB')


args = parser.parse_args()
dir_out=os.path.join('result',args.output_dir)
makedir(dir_out)
# result_mol_pic_out=os.path.join(dir_out,'result')
# makedir(result_mol_pic_out)
# result_file=os.path.join(args.result_dir,'result.csv')
# result_df=pd.read_csv(result_file)
# smiles=list(result_df.iloc[:,0])
# for i in range(len(smiles)):
#     s=smiles[i]
#     mol_png_file=os.path.join(result_mol_pic_out,f'{i}.png')
#     img = Draw.MolsToGridImage(mols=[Chem.MolFromSmiles(s) ],molsPerRow=1)
#     img.save(mol_png_file)





for iter in range(1,args.iter):
    iter_dir=os.path.join(dir_out,f'iter{iter}')
    makedir(iter_dir)
    fragment_file=os.path.join(args.result_dir,f'fragment_{iter}.csv')
    fragment_parent_file = os.path.join(args.result_dir, f'fragment_parent_{iter}.csv')
    df_data=pd.read_csv(fragment_file)
    df_parent=pd.read_csv(fragment_parent_file)
    fragment_list=list(df_data.columns)
    fragment_score_list=[]
    smiles_score_list=[]
    parent_list=[]
    for i in range(len(fragment_list)):
        score=df_data.iloc[0,i]
        score_list=list(df_data.iloc[1:,i])
        smiles_parent_list=list(df_parent.iloc[:,i])
        fragment_score_list.append(score)
        smiles_score_list.append(score_list)
        parent_list.append(smiles_parent_list)
    fragment_score_np=np.array(fragment_score_list)
    indx = list(fragment_score_np.argsort()[-args.num_fragment:])
    result_fragment=[]
    result_fragment_score=[]
    result_avg_score=[]
    for i in range(len(indx)) :
        id=indx[i]
        result_fragment.append(fragment_list[id])
        result_fragment_score.append(fragment_score_list[id])
        smiles_score_np=np.array(smiles_score_list[id])
        result_avg_score.append(np.nanmean(smiles_score_np))
        smiles_draw=[fragment_list[id]]+parent_list[id]
        # smiles_draw=[fragment_list[id]]
        scores_draw=[fragment_score_list[id]]+smiles_score_list[id]
        # scores_draw=[fragment_score_list[id]]
        smiles_draw_filted=[]
        scores_draw_filted=[]
        for j in range(len(scores_draw)):
            if not pd.isnull(smiles_draw[j]):
                smiles_draw_filted.append(smiles_draw[j])
                scores_draw_filted.append(scores_draw[j])
        png_fragment_and_parent_score=os.path.join(iter_dir,f'{i+1}.png')
        img = Draw.MolsToGridImage(mols=[Chem.MolFromSmiles(x) for x in smiles_draw_filted],legends=['%.1f'%x for x in scores_draw_filted], molsPerRow=len(scores_draw_filted))
        img.save(png_fragment_and_parent_score)
    result_df=pd.DataFrame({'fragment':result_fragment,'fragment_score':result_fragment_score,'avg_score':result_avg_score})
    result_file=os.path.join(dir_out,f'iter_{iter}.csv')
    result_df.to_csv(result_file,index=False)
    png_path = os.path.join(dir_out,f'iter_{iter}.png')
    img = Draw.MolsToGridImage(mols= [Chem.MolFromSmiles(x) for x in result_fragment],molsPerRow=10)
    img.save(png_path)



