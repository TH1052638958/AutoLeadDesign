import statistics

import openai
import jsonlines
import time
import pandas as pd
import numpy as np
import random
import json
from rdkit.Chem import Descriptors, QED
from rdkit.Chem import rdMolDescriptors as rdmd
# Set up your API key and model parameters



import openai
import jsonlines
import time
import subprocess
import re
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem import RDConfig
from rdkit.Chem import QED
import os
import sys
# from rdkit.Contrib.SA_Score import sascorer
#change to pre line
# sys.path.append(os.path.join(RDConfig.RDContribDir, 'SA_Score'))
# import sascorer
#change end
from tools.docking import calc_affinity
from tools.tools import get_fragment
import argparse
from openai import OpenAI
os.environ["OPENAI_BASE_URL"] = "https://key.wenwen-ai.com/v1"
os.environ["OPENAI_API_KEY"] = 'sk-I98mfNrsTBSVr4hfFa662506C3Fc4262B85f853dB45795B2'
sample=['C[C@@]1(O)CCc2ccccc21','N[C@@H](C)[C@@H]1COc2ccccc21','c1ccc2c(c1)NC(=O)CS2']
target_mol='C[C@@H]1C[C@H]2N[C@H]3Cc4cc(Cl)ccc4C[C@@H]3[C@@H]2N1'
# prompt = f'Generate a novel valid molecule SMILES which contains one fragment of [ {sample[0]} , {sample[1]},{sample[2]} ] at least and do not generate any English text.'
#prompt = f'Optimizate the molecule: {target_mol} to be more synthesizable and do not generate any English text'
prompt=f'I have three molecular fragments represented in SMILES notation: [SMILES1], [SMILES2], and [SMILES3]. Please generate a new molecule in SMILES format that incorporates at least one of these fragments. The generated molecule should be realistic and chemically plausible, with a structure that can potentially be synthesized. Aim for a diverse and creative structure while ensuring that it contains at least one of the given fragments.'
client = OpenAI()
completion = client.chat.completions.create(
                model='gpt-4-turbo',
                messages=[
                    {"role": "user", "content": prompt}],
                n=1,
                max_tokens=60,
                temperature=0.6,
                stop="!",
                user="user"
            )
new_mol=json.loads(completion.model_dump_json())["choices"][0]["message"]["content"]
print(new_mol)
score=calc_affinity(new_mol,dir_out='smina',name_protein='6igx')
print(score)