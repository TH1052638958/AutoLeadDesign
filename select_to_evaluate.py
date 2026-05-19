import rdkit
from rdkit import Chem

m=Chem.MolFromSmiles('CC1=C(NC(=O)C2=CC=CC(C)=C2C2=NC(N)=CC=C2)C=CC=C1NCC1=CC=C2C(=C1)C1=NC(N)=CC=C1C2CC=N')
ri=m.GetRingInfo()

for ring in ri.AtomRings():
    print(ring)
    if len(ring) == 1:
        print('1')