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



mol=Chem.MolFromSmiles('n1c(=O)ncc2cc(F)cnc21')
drawer = rdMolDraw2D.MolDraw2DCairo(300, 300)  # 可以自行调整图像尺寸
drawer.drawOptions().clearBackground = False
drawer.drawOptions().addStereoAnnotation = False
drawer.DrawMolecule(mol)
drawer.FinishDrawing()
png = drawer.GetDrawingText()
bio = io.BytesIO(png)
img = Image.open(bio)
