import pandas as pd
import os
import argparse

parser = argparse.ArgumentParser(description='LMLF')
parser.add_argument('--reference', type=str, default='random')
parser.add_argument('--protein-name', type=str, default='1')
parser.add_argument('--method', type=str, default='ALD')
args = parser.parse_args()
result_file=os.path.join(args.method+'_'+args.reference+'_'+args.protein_name,'result.csv')
data_df=pd.read_csv(result_file)
score=list(data_df.iloc[:,1])
score.sort(reverse=True)
top1=score[0]
top10=score[0:10]
top100=score[0:100]
print('Top1:',top1)
print('Top1:',sum(top10)/len(top10))
print('Top1:',sum(top100)/len(top100))
