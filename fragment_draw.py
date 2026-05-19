import os.path

from rdkit import Chem
from rdkit import rdBase
from rdkit.Chem import Draw
from rdkit.Chem.Draw import rdMolDraw2D
from rdkit.Chem import rdDepictor
rdDepictor.SetPreferCoordGen(True)
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem import rdFMCS
from rdkit import DataStructs
import io
from PIL import Image
from tools.tools import makedir

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

import  pandas as pd
from tools.tools import get_fragment
# exp_path='fragment_PRMT5'
# data=pd.read_csv(os.path.join(exp_path,'init_score.csv'))
# smiles_list=list(data.iloc[:,0])
smiles_list=['CC1CN(CCN1C2=NC(=O)N(C3=NC(=C(C=C32)F)C4=C(C=CC=C4F)O)C5=C(C=CN=C5C(C)C)C)C(=O)C=C']
exp_path='fig'
for i in range(len(smiles_list)):
    save_dir=os.path.join(exp_path,str(i))
    makedir(save_dir)
    smiles=smiles_list[i]
    fragment=get_fragment(smiles,'BRICS')
    png=get_png(smiles)
    png.save(os.path.join(save_dir,'smiles.png'))
    for j in range(len(fragment)):
        save_file=os.path.join(save_dir,str(j+1)+'.png')
        f=list(fragment)[j]
        png=get_png(f)
        png.save(save_file)

