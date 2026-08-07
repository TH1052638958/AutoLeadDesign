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
from tools.tools import draw_smiles
from tools.tools import makedir
from rdkit.Chem import Draw

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


dir_name='fragment_algorithm'
makedir(dir_name)

smiles='O=C(COc1ccc(CO)cc1)N[C@@H]1[C@@H]2Cc3cc(Cl)ccc3[C@@H]21'
smiles_list=['CCNC(=O)c1ccc(NC[C@@]2(O)CCc3ccccc32)nc1','CN(Cc1nc2ccccc2s1)C(=O)[C@H]1Cc2ccccc2O1','O=C(N[C@H]1CCCNC1=O)c1cc(CCc2ccccc2)ccc1O'
             ,'Cc1ccc([C@H]2C[C@H]2C(=O)NCc2ccc3nc(O)[nH]c3c2)cc1','O=C(COc1ccc(CO)cc1)N[C@@H]1[C@@H]2Cc3cc(Cl)ccc3[C@@H]21'
             ,'O=C(CC/C(O)=N/[C@@H]1CCCN=C1O)N1CCSc2ccccc21','NC(=O)C1=NO[C@@H](CNC(=O)N2CCC(c3ccc(O)cc3)CC2)C1']
algorithm_list=['BRICS','RECAP']
for i in range(len(smiles_list)):
    exp_path=os.path.join(dir_name,str(i))
    makedir(exp_path)
    for j in algorithm_list:
        fragment=get_fragment(smiles_list[i],j)
        png_path=os.path.join(exp_path,j+'.png')
        img=Draw.MolsToGridImage(mols=[Chem.MolFromSmiles(smiles_list[i])]+[Chem.MolFromSmiles(x) for x in fragment])
        img.save(png_path)







