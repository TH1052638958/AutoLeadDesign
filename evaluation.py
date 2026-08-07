from rdkit import Chem
from rdkit.Chem import AllChem, DataStructs
import numpy as np
import pandas as pd
def calculate_internal_diversity(smiles_list):
    diversity_matrix = calculate_tanimoto_diversity(smiles_list)
    n = len(smiles_list)
    total_diversity = np.sum(diversity_matrix) / (n * (n - 1) / 2)
    return total_diversity

def calculate_tanimoto_diversity(smiles_list):
    fps = [AllChem.GetMorganFingerprintAsBitVect(Chem.MolFromSmiles(s), 2) for s in smiles_list]
    n = len(fps)
    diversity_matrix = np.zeros((n, n))

    for i in range(n):
        for j in range(i + 1, n):
            sim = DataStructs.TanimotoSimilarity(fps[i], fps[j])
            diversity_matrix[i, j] = 1 - sim
            diversity_matrix[j, i] = 1 - sim

    return diversity_matrix
rear_df=pd.read_csv('result_min.csv')
fragment_df=pd.read_csv('fragment_Ampc/result.csv')
smiles_rear=list(rear_df.iloc[:,0])
smiles_fragment=list(fragment_df.iloc[:,0])
print(smiles_rear[0])
# diversity_rear=calculate_tanimoto_diversity(smiles_rear)
# diversity_fragment=calculate_tanimoto_diversity(smiles_fragment)
diversity_rear=calculate_internal_diversity(smiles_rear)
diversity_fragment=calculate_internal_diversity(smiles_fragment)
print('rear:',diversity_rear)
print('fragment:',diversity_fragment)



