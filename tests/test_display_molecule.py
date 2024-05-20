import unittest
from unittest.mock import patch, MagicMock
import sys

# Adjust the path to include the directory containing your Chem_pack module
sys.path.insert(0, "./src")

from project_ppchem_tools_kit.display_molecule import display_molecule

class TestDisplayMolecule(unittest.TestCase):
    def setUp(self):
        # Patch the Chem module
        self.patcher_chem = patch('project_ppchem_tools_kit.display_molecule.Chem.MolFromSmiles', autospec=True)
        self.mock_mol_from_smiles = self.patcher_chem.start()

        # Patch the Draw module
        self.patcher_draw = patch('project_ppchem_tools_kit.display_molecule.Draw.MolToImage', autospec=True)
        self.mock_mol_to_image = self.patcher_draw.start()

        # Patch the Tkinter components
        self.patcher_tk = patch('project_ppchem_tools_kit.display_molecule.tk.Toplevel', autospec=True)
        self.mock_toplevel = self.patcher_tk.start()

        # Mock for FigureCanvasTkAgg
        self.patcher_canvas = patch('project_ppchem_tools_kit.display_molecule.FigureCanvasTkAgg', autospec=True)
        self.mock_canvas = self.patcher_canvas.start()

        # Mock for NavigationToolbar2Tk
        self.patcher_toolbar = patch('project_ppchem_tools_kit.display_molecule.NavigationToolbar2Tk', autospec=True)
        self.mock_toolbar = self.patcher_toolbar.start()

        # Create a mock figure to be returned by FigureCanvasTkAgg
        self.mock_figure = MagicMock()
        self.mock_canvas_instance = MagicMock()
        self.mock_canvas.return_value = self.mock_canvas_instance
        self.mock_canvas_instance.figure = self.mock_figure

    def tearDown(self):
        self.patcher_chem.stop()
        self.patcher_draw.stop()
        self.patcher_tk.stop()
        self.patcher_canvas.stop()
        self.patcher_toolbar.stop()

    def test_display_molecule_valid_smiles(self):
        # Mock entry_input to return a valid SMILES string
        mock_entry_input = MagicMock()
        mock_entry_input.get.return_value = 'C1=CC=CC=C1'  # Benzene SMILES

        # Mock MolFromSmiles to return a valid molecule object
        mock_mol = MagicMock()
        self.mock_mol_from_smiles.return_value = mock_mol

        # Mock MolToImage to return a dummy image
        mock_image = MagicMock()
        self.mock_mol_to_image.return_value = mock_image

        # Call the function
        display_molecule(mock_entry_input)

        # Assertions to check if the functions were called correctly
        self.mock_mol_from_smiles.assert_called_once_with('C1=CC=CC=C1')
        self.mock_mol_to_image.assert_called_once_with(mock_mol)
        self.mock_toplevel.assert_called_once()

    def test_display_molecule_invalid_smiles(self):
        # Mock entry_input to return an invalid SMILES string
        mock_entry_input = MagicMock()
        mock_entry_input.get.return_value = 'invalid_smiles'

        # Mock MolFromSmiles to return None for an invalid SMILES string
        self.mock_mol_from_smiles.return_value = None

        # Call the function and capture print output
        with patch('builtins.print') as mocked_print:
            display_molecule(mock_entry_input)
            mocked_print.assert_called_once_with("Erreur : Impossible de générer une structure moléculaire à partir du SMILES fourni.")

if __name__ == '__main__':
    unittest.main()
