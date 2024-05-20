
import pubchempy as pcp


def name_to_smiles(molecule_name):
    """Fetches the SMILES representation of a molecule from its name from PubChem.

    Args:
        molecule_name (str): The name of the molecule.

    Returns:
        str: The canonical SMILES string if found, None otherwise.
    
    Raises:
        PubChemHTTPError: If any issue is encountered with the PubChem API request.
    """
    try:
        compound = pcp.get_compounds(molecule_name, 'name')
        if compound:
            return compound[0].canonical_smiles
        else:
            return None
    except pcp.PubChemHTTPError as e:
        print("Error occurred while fetching data from PubChem:", e)
        return None
    

