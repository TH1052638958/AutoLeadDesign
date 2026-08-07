import  os
from rdkit.Chem import Descriptors, QED
from rdkit.Chem import rdMolDescriptors as rdmd
from rdkit import Chem
# from rdkit.Contrib.SA_Score import sascorer
from rdkit.Chem import RDConfig
import os
import sys
sys.path.append(os.path.join(RDConfig.RDContribDir, 'SA_Score'))
# now you can import sascore!
import sascorer
from rdkit.Chem.Draw import rdMolDraw2D
import io
from PIL import Image
def makedir(path):
    try:
        os.makedirs(path)
    except OSError:
        pass

"""
Lipinski
"""
from rdkit.Chem import Descriptors, Crippen, Lipinski
from copy import deepcopy
import numpy as np
import re
def extract_SMILES(text):
    result_split=text.split()
    smiles=None
    for s in result_split[::-1]:
        if len(s)>=7:
            try:
                mol=Chem.MolFromSmiles(s)
                smiles=s
                if mol==None:
                    continue
                break
            except Exception as e:
                continue

    return smiles


def obey_lipinski(smile):
    mol = Chem.MolFromSmiles(smile)
    Chem.SanitizeMol(mol)
    rule_1 = Descriptors.ExactMolWt(mol) < 500
    rule_2 = Lipinski.NumHDonors(mol) <= 5
    rule_3 = Lipinski.NumHAcceptors(mol) <= 10
    rule_4 = (logp:=Crippen.MolLogP(mol)>=-2) & (logp<=5)
    rule_5 = Chem.rdMolDescriptors.CalcNumRotatableBonds(mol) <= 10
    return np.sum([int(a) for a in [rule_1, rule_2, rule_3, rule_4, rule_5]])






def calculate_radscore(mol):
    # molecule = Chem.MolFromSmiles(mol)
    # sa_score = rdMolDescriptors.SyntheticAccessibility(molecule)
    try:
        m = Chem.MolFromSmiles(mol)
        sa_score = sascorer.calculateScore(m)
        return sa_score
    except:
        return None
def calculate_qed_score(smiles):
    molecule = Chem.MolFromSmiles(smiles)
    if molecule is not None:
        qed_score = QED.qed(molecule)
        return qed_score
    else:
        return None


def Sdf2SMILES(base_dir):
    result=[]
    file_list=[os.path.join(base_dir, file) for file in os.listdir(base_dir)]
    for f in file_list:
        mol=Chem.SDMolSupplier(f)
        SMI=Chem.MolToSmiles(mol[0])
        result.append(SMI)
    return result




def get_fragment(smiles,method,process=True):
    if method=='BRICS':
        return extract_fragment_by_BRICS(smiles,process)
    elif method=='RECAP':
        return extract_fragment_by_RECAP(smiles,process)

def extract_fragment_by_BRICS(smiles,process):
    from rdkit.Chem import BRICS
    from rdkit import Chem
    aspirin = Chem.MolFromSmiles(smiles)
    fragments = BRICS.BRICSDecompose(aspirin, allNodes=None, minFragmentSize=2,
                                     onlyUseReactions=None, silent=True, keepNonLeafNodes=False, singlePass=False,
                                     returnMols=True)
    fragment_mol=fragments
    fragment_s=[]
    for i in fragment_mol:

        fragment_s.append(Chem.MolToSmiles(i,kekuleSmiles=True))
    fragments=fragment_s
    if process:

        #remove [5*] and ([5*])
        result=[]
        for f in fragments:
            for i in range(1,21):
                f=f.replace(f'([{i}*])','')
                f=f.replace(f'[{i}*]','')
            result.append(f)

    else:
        result=fragments
    return result

def extract_fragment_by_RECAP(smiles,process):
    from rdkit.Chem import Recap
    m = Chem.MolFromSmiles(smiles)
    hierarch = Recap.RecapDecompose(m)
    result=[Chem.MolToSmiles(x.mol) for x in hierarch.GetLeaves().values()]
    #remove * and (*)
    if process :

        result_filted=[]
        for i in result:
            i=i.replace('(*)','')
            i=i.replace('*','')
            result_filted.append(i)
        result=result_filted
    return result


def get_png(smiles):
    mol=Chem.MolFromSmiles(smiles)
    drawer = rdMolDraw2D.MolDraw2DCairo(300, 300)  # 可以自行调整图像尺寸
    drawer.drawOptions().clearBackground = False
    drawer.drawOptions().addStereoAnnotation = False
    drawer.DrawMolecule(mol)
    drawer.FinishDrawing()
    png = drawer.GetDrawingText()
    bio = io.BytesIO(png)
    img = Image.open(bio)
    # img.save('result.png')
    return img

def draw_smiles(smiles,file_path):
    png = get_png(smiles)
    png.save(file_path)

def get_fingerpoint_by_MACCS(smiles_list):
    from rdkit.Chem import AllChem
    fp_list=[]
    for smiles in smiles_list:

        try:
            mol = Chem.MolFromSmiles(smiles)
            fp=AllChem.GetMACCSKeysFingerprint(mol)
            fp=fp.ToBitString()
        except Exception as e:
            print(f'cont gen MACCS for {smiles}',e)
            fp=None
        if fp !=None:
            fp_list.append(fp)
    return fp_list

def get_fingerpoint_by_Topological(smiles_list):
    from rdkit.Chem import AllChem
    from rdkit.Chem.Fingerprints import FingerprintMols
    fp_list=[]
    for smiles in smiles_list:
        mol=Chem.MolFromSmiles(smiles)
        try:
            fp=FingerprintMols.FingerprintMol(mol)
            fp=fp.ToBitString()
        except Exception as e:
            print(f'cont gen MACCS for {smiles}',e)
            fp=None
        fp_list.append(fp)
    return fp_list

def get_fingerpoint_by_Morgan(smiles_list):
    from rdkit.Chem import AllChem
    fp_list=[]
    for smiles in smiles_list:
        mol=Chem.MolFromSmiles(smiles)
        try:
            fp=AllChem.GetMorganFingerprintAsBitVect(mol, 2, 2048)
            fp=fp.ToBitString()
        except Exception as e:
            print(f'cont gen MACCS for {smiles}',e)
            fp=None
        fp_list.append(fp)
    return fp_list

def get_fingerpoint_by_Avalon(smiles_list):
    from rdkit.Avalon import pyAvalonTools
    fp_list=[]
    for smiles in smiles_list:
        mol=Chem.MolFromSmiles(smiles)
        try:
            fp=pyAvalonTools.GetAvalonFP(mol,nBits=128)
            fp=fp.ToBitString()
        except Exception as e:
            print(f'cont gen MACCS for {smiles}',e)
            fp=None
        fp_list.append(fp)
    return fp_list


def get_fingerpoint(smiles_list,method):
    if method=='MACCS':
        return get_fingerpoint_by_MACCS(smiles_list)
    elif method=='Topological':
        return get_fingerpoint_by_Topological(smiles_list)
    elif method=='Morgan':
        return get_fingerpoint_by_Morgan(smiles_list)
    elif method=='Avalon':
        return get_fingerpoint_by_Avalon(smiles_list)
    else:
        raise ValueError('Method Error')