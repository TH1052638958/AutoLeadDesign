from rdkit import Chem
from rdkit.Chem import MACCSkeys
from rdkit.Chem.AtomPairs import Torsions
def get_fingerpoint(smiles_list):
    result=[]
    for smi in smiles_list:
        mol= Chem.MolFromSmiles(smi)
        fp = Torsions.GetTopologicalTorsionFingerprintAsIntVect(mol)

        result.append(fp.ToBinary())




