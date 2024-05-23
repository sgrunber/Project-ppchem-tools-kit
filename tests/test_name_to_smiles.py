import unittest
from unittest.mock import patch

import sys

from project_ppchem_tools_kit.name_to_smiles import name_to_smiles

#2.NAMETOSMILES
class TestChemFunctions(unittest.TestCase):
    @patch('project_ppchem_tools_kit.name_to_smiles.pcp.get_compounds')
    def test_name_to_smiles(self, mock_get_compounds):
        # Setup
        mock_get_compounds.return_value = [type('test', (object,), {"canonical_smiles": "C"})()]  # Simulated response
        # Execution & Assertion
        self.assertEqual(name_to_smiles("water"), "C")

        # Test for a None response
        mock_get_compounds.return_value = []
        self.assertIsNone(name_to_smiles("unknown"))
if __name__ == '__main__':
    unittest.main()