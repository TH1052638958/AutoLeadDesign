import argparse
import pandas as pd
import os
import json
import jsonlines
parser = argparse.ArgumentParser(description='LMLF')
parser.add_argument('--exp-path', type=str, default='6IGX')
parser.add_argument('--init-path', type=str, default='6igx')
parser.add_argument('--save', type=int, default=1)

args = parser.parse_args()
init_csv=os.path.join(args.init_path,'init_score.csv')
init_jsonl=os.path.join(args.init_path,'init.jsonl')

data_df=pd.read_csv(init_csv)
smiles=list(data_df.iloc[:,0])
scores=list(data_df.iloc[:,1])
threshold=int(scores[-1]*0.9)
if args.save==0:
    print(threshold)
else:
    smiles=smiles[:-1]
    scores=scores[:-1]
    result_csv=os.path.join(args.exp_path,'init_score.csv')
    result1_csv = os.path.join(args.exp_path, 'init.csv')
    result_jsonl=os.path.join(args.exp_path,'init.jsonl')
    pd.DataFrame({'smile':smiles,'score':scores}).to_csv(result_csv,index=False)
    labels=[]
    with jsonlines.open(init_jsonl) as reader:
        for line in reader:
            if "\n" not in line:
                labels.append(line['label'])
    json_result=[]
    for i in range(len(smiles)):
        json_result.append({'smiles': smiles[i], 'label': labels[i],'score':str(scores[i])})
    with jsonlines.open(result_jsonl, mode='a') as writer:
        # writer.write("\\n")
        for molecule in json_result:
            writer.write(molecule)
            writer.write('\n')
    pd.DataFrame({'smile': smiles}).to_csv(result1_csv, index=False)

    print(threshold)