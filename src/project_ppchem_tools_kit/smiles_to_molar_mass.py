from rdkit import Chem
from rdkit.Chem import Descriptors

def smiles_to_molar_mass(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is not None:
        return Descriptors.ExactMolWt(mol)
    else:
        return None
    