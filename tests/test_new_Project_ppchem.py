import sys
sys.path.insert(0, "./src")

#importing functions from src/Chem_pack 
from Chem_pack.smiles_to_molar_mass import smiles_to_molar_mass
from Chem_pack.name_to_smiles import name_to_smiles
from Chem_pack.display_molecule import display_molecule
from Chem_pack.make_graph import make_graph
from Chem_pack.linear_regression import linear_regression
from Chem_pack.error_calculation_interface import error_calculation_interface
from Chem_pack.process_input import process_input
from Chem_pack.copy_text import copy_text
from Chem_pack.select_all import select_all

import unittest
from unittest.mock import patch, MagicMock
import tkinter as tk

#importing functions from the project file
'''import Project_ppchem
from Project_ppchem import relative_to_assets, ASSETS_PATH, name_to_smiles, smiles_to_molar_mass, linear_regression, make_graph, process_input, display_result, browse_excel_file, select_all, copy_text, welcome_message, clear_input, on_radio_select, create_radio_button, bind_enter
'''

from pathlib import Path
import os
OUTPUT_PATH = Path(os.getcwd())
ASSETS_PATH = OUTPUT_PATH / Path("./assets/frame0")

def relative_to_assets(path: str) -> Path:
    return ASSETS_PATH / Path(path)

#1.RELATIVETOASSETS
class TestRelativePath(unittest.TestCase):
    def test_relative_to_assets(self):
        # Test the path computation is correct
        test_path = "subfolder/file.txt"
        expected = ASSETS_PATH / "subfolder/file.txt"
        result = relative_to_assets(test_path)
        self.assertEqual(result, expected)
if __name__ == '__main__':
    unittest.main()

'''if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTest(unittest.makeSuite(TestRelativePath))
    runner = unittest.TextTestRunner()
    runner.run(suite)
'''

#2.NAMETOSMILES
class TestChemFunctions(unittest.TestCase):
    @patch('Chem_pack.name_to_smiles.pcp.get_compounds')
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

'''if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTest(unittest.makeSuite(TestChemFunctions))
    runner = unittest.TextTestRunner()
    runner.run(suite)
'''

#3.SMILESTOMOLARMASS
class TestMolecularFunctions(unittest.TestCase):
    @patch('Chem_pack.smiles_to_molar_mass.Chem.MolFromSmiles')
    @patch('Chem_pack.smiles_to_molar_mass.Descriptors.ExactMolWt')
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

'''if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTest(unittest.makeSuite(TestMolecularFunctions))
    runner = unittest.TextTestRunner()
    runner.run(suite)
'''


#4,5 LINEARREG & MAKEGRAPH
import unittest
from unittest.mock import Mock, patch
from Chem_pack.linear_regression import linear_regression
from Chem_pack.make_graph import make_graph
import pandas as pd 

class TestPlottingFunctions(unittest.TestCase):
    @patch('Chem_pack.linear_regression.pd.read_excel')
    @patch('Chem_pack.linear_regression.display_graph')
    def test_linear_regression(self, mock_display_graph, mock_read_excel):
        # Mock the return value of pd.read_excel
        mock_read_excel.return_value = pd.DataFrame({
            'X': [1, 2, 3, 4, 5],
            'Y': [2, 4, 6, 8, 10]
        })
        
        # Create a mock for entry_input
        mock_entry_input = Mock()
        mock_entry_input.get.return_value = 'mock_path.xlsx'
        
        # Call the function with the mock object
        linear_regression(mock_entry_input, "X", "Y", "Test")

        # Check if display_graph was called once
        mock_display_graph.assert_called_once()

    @patch('Chem_pack.make_graph.pd.read_excel')
    @patch('Chem_pack.make_graph.display_graph')
    def test_make_graph(self, mock_display_graph, mock_read_excel):
        # Mock the return value of pd.read_excel
        mock_read_excel.return_value = pd.DataFrame({
            'X': [1, 2, 3, 4, 5],
            'Y': [2, 4, 6, 8, 10]
        })
        
        # Create a mock for entry_input
        mock_entry_input = Mock()
        mock_entry_input.get.return_value = 'mock_path.xlsx'
        
        # Call the function with the mock object
        make_graph(mock_entry_input, "X", "Y", "Test")

        # Check if display_graph was called once
        mock_display_graph.assert_called_once()

if __name__ == '__main__':
    unittest.main()



'''
#6 PROCESS INPUT
class TestProcessInput(unittest.TestCase):
    def setUp(self):
        # Patching the Entry widget
        self.patcher_entry = patch('tkinter.Entry', autospec=True)
        self.mock_entry = self.patcher_entry.start()
        self.mock_entry_instance = MagicMock()
        self.mock_entry.return_value = self.mock_entry_instance
        self.mock_entry_instance.get.return_value = ''

        # Mock the messagebox to prevent actual GUI dialog boxes during the test
        self.patcher_msgbox = patch('tkinter.messagebox', autospec=True)
        self.mock_messagebox = self.patcher_msgbox.start()

    def tearDown(self):
        self.patcher_entry.stop()
        self.patcher_msgbox.stop()

    def test_empty_input(self):
        # Ensure entry_input is mocked
        with patch('Chem_pack.process_input.entry_input', self.mock_entry_instance):
            process_input()
            self.mock_messagebox.showerror.assert_called_once_with("Error", "Please enter a molecule name, SMILES code, or file path.")
if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTest(unittest.makeSuite(TestProcessInput))
    runner = unittest.TextTestRunner()
    runner.run(suite)



#9. SELECTALL
import tkinter as tk
class TestSelectAll(unittest.TestCase):
    def test_select_all(self):
        """Test select_all functionality."""
        mock_event = MagicMock()
        mock_widget = MagicMock()
        mock_event.widget = mock_widget
        select_all(mock_event)
        mock_widget.tag_add.assert_called_with("sel", "1.0", "end")
if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTest(unittest.makeSuite(TestSelectAll))
    runner = unittest.TextTestRunner()
    runner.run(suite)

#10. COPYTEXT
class TestCopyText(unittest.TestCase):
    def test_copy_text(self):
        """Test copy_text triggers the copy event."""
        mock_event = MagicMock()
        mock_widget = MagicMock()
        mock_event.widget = mock_widget
        copy_text(mock_event)
        mock_widget.event_generate.assert_called_with("<<Copy>>")
if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTest(unittest.makeSuite(TestCopyText))
    runner = unittest.TextTestRunner()
    runner.run(suite)

'''