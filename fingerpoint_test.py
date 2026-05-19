


# import os
# from rdkit import Chem
# from rdkit.Chem import Draw
# from rdkit.Chem import RDConfig
#
# from rdkit.Chem import FragmentCatalog
#
# fName = os.path.join(RDConfig.RDDataDir, 'FunctionalGroups.txt')
#
# fparams = FragmentCatalog.FragCatParams(1, 6, fName)
# print(fparams.GetNumFuncGroups())
# for i in range(fparams.GetNumFuncGroups()):
#     m = fparams.GetFuncGroup(i)
#
#     print(Chem.MolToSmiles(m))



#MACCS
# from rdkit import Chem
#
# mol = Chem.MolFromSmiles('CC(C)C1=C(C(=C(N1CC[C@H](C[C@H](CC(=O)O)O)O)C2=CC=C(C=C2)F)C3=CC=CC=C3)C(=O)NC4=CC=CC=C4')
# from rdkit.Chem import MACCSkeys
#
# fp = MACCSkeys.GenMACCSKeys(mol)
# print(type(fp))  # <class 'rdkit.DataStructs.cDataStructs.ExplicitBitVect'>
# for i in range(len(fp)):
#     print(fp[i], end='')
# print('/n')
# fp_bits = tuple(fp.GetOnBits())
# print(fp_bits)



#药效团

# from rdkit import Chem
# from rdkit.Chem import ChemicalFeatures
# from rdkit import RDConfig
# import os
# fdefName = os.path.join(RDConfig.RDDataDir,'BaseFeatures.fdef')
# factory = ChemicalFeatures.BuildFeatureFactory(fdefName)
# from rdkit.Chem.Pharm2D.SigFactory import SigFactory
# featFactory = ChemicalFeatures.BuildFeatureFactory(fdefName)
# sigFactory = SigFactory(featFactory, minPointCount=2, maxPointCount=3)
# sigFactory.SetBins([(0,2), (2,5), (5,8)])
# sigFactory.Init()
# print(sigFactory.GetSigSize())
# from rdkit.Chem.Pharm2D import Generate
# mol = Chem.MolFromSmiles('OCC(=O)CCCN')
# fp = Generate.Gen2DFingerprint(mol,sigFactory)
# for bit in range(sigFactory.GetSigSize()):
#     if fp.GetBit(bit):
#         print(f"Bit {bit}: {sigFactory.GetBitDescription(bit)}")
# print(len(fp))
# print(fp.GetNumOnBits())
# print(sigFactory.GetBitDescription(100))
# print(sigFactory.GetFeatFamilies())
from rdkit.Chem import AllChem


#断键1
# from rdkit.Chem import BRICS
# from rdkit import Chem
# aspirin= Chem.MolFromSmiles('CC(=O)OC1=CC=CC=C1C(O)=O')
# fragments=BRICS.BRICSDecompose(aspirin,allNodes=None, minFragmentSize=2,
# onlyUseReactions=None, silent=True, keepNonLeafNodes=False, singlePass=False, returnMols=False)
# print(fragments)


#断键
# from rdkit import Chem
# from rdkit.Chem import BRICS
# import numpy as np
# from rdkit.Chem import MolToSmiles as mol_to_smiles
# def fragment_recursive(mol, frags):
#     try:
#         bonds = list(BRICS.FindBRICSBonds(mol))
#         if len(bonds) == 0:
#             frags.append(mol_to_smiles(mol))
#             return frags
#         idxs, labs = list(zip(*bonds))
#         bond_idxs = []
#         for a1, a2 in idxs:
#             bond = mol.GetBondBetweenAtoms(a1, a2)
#             bond_idxs.append(bond.GetIdx())
#         order = np.argsort(bond_idxs).tolist()
#         bond_idxs = [bond_idxs[i] for i in order]
#         broken = Chem.FragmentOnBonds(mol,
#                                       bondIndices=[bond_idxs[0]],
#                                       dummyLabels=[(0, 0)])
#         head, tail = Chem.GetMolFrags(broken, asMols=True)
#         #print(mol_to_smiles(head), mol_to_smiles(tail))
#         frags.append(mol_to_smiles(head))
#         return fragment_recursive(tail, frags)
#     except Exception as e:
#         print (e)
#         pass
#
# aspirin= Chem.MolFromSmiles('CC(=O)OC1=CC=CC=C1C(O)=O')
# fragments=fragment_recursive(aspirin, [])
# print (fragments)

from rdkit import Chem
from rdkit.Chem import rdmolops



#药效团 success
import os
from rdkit import Geometry
from rdkit import RDConfig
from rdkit.Chem import AllChem
from rdkit.Chem import ChemicalFeatures
from rdkit import Chem
from rdkit.Chem import Draw
import matplotlib.pyplot as plt
from rdkit.Chem.Pharm3D import Pharmacophore
def extract_fragment(smiles, atom_ids):
    # 创建分子对象
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None

    # 获取原子索引
    atom_ids = [int(i) for i in atom_ids]

    # 创建子结构
    fragment = Chem.PathToSubmol(mol, atom_ids)

    # 返回片段的 SMILES
    return Chem.MolToSmiles(fragment)

FEAT = os.path.join(RDConfig.RDDataDir, "BaseFeatures.fdef")
featfact = ChemicalFeatures.BuildFeatureFactory(FEAT)
smile='CC(=O)OC1=CC=CC=C1C(O)=O'

mol = Chem.MolFromSmiles(smile)
for atom in mol.GetAtoms():
    atom.SetAtomMapNum(atom.GetIdx())

AllChem.EmbedMolecule(mol)
feats = featfact.GetFeaturesForMol(mol)
for feat in feats:
    print(feat.GetFamily())
    # pos = feat.GetPos()
    # print(pos.x, pos.y, pos.z)
    print(feat.GetAtomIds())
    img = Draw.MolToImage(mol, size=(300, 300), kekulize=True,highlightAtoms=list(feat.GetAtomIds()) ,wedgeBonds=True)
    plt.imshow(img)
    plt.axis('off')
    plt.show()
    print(extract_fragment(smile,feat.GetAtomIds()))




# from rdkit import Chem
# from rdkit.Chem import Draw
# import matplotlib.pyplot as plt


# def draw_molecule_with_atom_numbers(smiles):
#     # 将 SMILES 转换为分子对象
#     mol = Chem.MolFromSmiles(smiles)
#
#     # 生成二维坐标
#     # Chem.rdDepictor.Compute2DCoords(mol)
#
#     # 创建原子编号的标签
#     atom_labels = [str(atom.GetIdx()) for atom in mol.GetAtoms()]
#
#     # 绘制分子结构
#     img = Draw.MolToImage(mol, size=(300, 300), kekulize=True, wedgeBonds=True, atomLabel=atom_labels)
#
#     # 显示图像
#     plt.imshow(img)
#     plt.axis('off')  # 关闭坐标轴
#     plt.show()
#
#
# # 示例 SMILES
# smiles_input = "CCO"  # 乙醇的 SMILES
# draw_molecule_with_atom_numbers(smiles_input)

