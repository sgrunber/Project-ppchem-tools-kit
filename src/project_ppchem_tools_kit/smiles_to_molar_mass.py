from rdkit import Chem
from rdkit.Chem import Descriptors

def smiles_to_molar_mass(smiles):
    """Calculates the molar mass of a molecule from a given SMILES using RDKit.

    Args:
        smiles (str): The SMILES string of the molecule.

    Returns:
        float: The molar mass of the molecules (in grams per mole) if the SMILES string is valid, None otherwise.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is not None:
        return Descriptors.ExactMolWt(mol)
    else:
        return None
    