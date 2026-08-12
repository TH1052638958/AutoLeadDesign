import os
import argparse
import pandas as pd

parser = argparse.ArgumentParser()
parser.add_argument("--exp-path", type=str, required=True)
parser.add_argument("--result", type=str, default="init_score.csv")
parser.add_argument("--rate", type=float, default=0.7)
args = parser.parse_args()

result_path = os.path.join(args.exp_path, args.result)
data = pd.read_csv(result_path)

max_score = data["dcoking_scores"].max()
threshold = max_score * args.rate

print(threshold)