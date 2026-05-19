from tools.tools import extract_fragment_by_BRICS ,extract_fragment_by_RECAP


smiles='CC1=C(NC(=O)C2=CC=CC(C)=C2C2=NC(N)=CC=C2)C=CC=C1NCC1=CC=C2C(=C1)C1=NC(N)=CC=C1C2CC=N'
fragment=extract_fragment_by_BRICS(smiles,process=True)
print(fragment)
# smiles='CCNC(=O)c1ccc(NC[C@@]2(O)CCc3ccccc32)nc1'
# fragment=extract_fragment_by_BRICS(smiles)
# print(fragment)
