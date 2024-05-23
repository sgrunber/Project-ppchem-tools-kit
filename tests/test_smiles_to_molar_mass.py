import sys

from project_ppchem_tools_kit.smiles_to_molar_mass import smiles_to_molar_mass

import unittest
from unittest.mock import patch, MagicMock


#3.SMILESTOMOLARMASS
class TestMolecularFunctions(unittest.TestCase):
    @patch('project_ppchem_tools_kit.smiles_to_molar_mass.Chem.MolFromSmiles')
    @patch('project_ppchem_tools_kit.smiles_to_molar_mass.Descriptors.ExactMolWt')
    def test_smiles_to_molar_mass(self, mock_ExactMolWt, mock_MolFromSmiles):
        # Setup
        mock_mol = MagicMock()
        mock_MolFromSmiles.return_value = mock_mol
        mock_ExactMolWt.return_value = 18.015

        # Execution & Assertion
        self.assertEqual(smiles_to_molar_mass("H2O"), 18.015)

        # Test for a None SMILES
        mock_MolFromSmiles.return_value = None
        self.assertIsNone(smiles_to_molar_mass("invalid"))
if __name__ == '__main__':
    unittest.main()