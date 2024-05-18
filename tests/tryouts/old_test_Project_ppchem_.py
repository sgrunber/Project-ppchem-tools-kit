import unittest
from unittest.mock import patch, MagicMock
from ToolsKit.Project_ppchem import relative_to_assets, ASSETS_PATH


#RELATIVETOASSETS
class TestRelativePath(unittest.TestCase):
    def test_relative_to_assets(self):
        # Test the path computation is correct
        test_path = "subfolder/file.txt"
        expected = ASSETS_PATH / "subfolder/file.txt"
        result = relative_to_assets(test_path)
        self.assertEqual(result, expected)

if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTest(unittest.makeSuite(TestRelativePath))
    runner = unittest.TextTestRunner()
    runner.run(suite)


#NAMETOSMILES
from ToolsKit.Project_ppchem import name_to_smiles

class TestChemFunctions(unittest.TestCase):
    @patch('Project_ppchem.pcp.get_compounds')
    def test_name_to_smiles(self, mock_get_compounds):
        # Setup
        mock_get_compounds.return_value = [type('test', (object,), {"canonical_smiles": "C"})()]  # Simulated response
        # Execution & Assertion
        self.assertEqual(name_to_smiles("water"), "C")

        # Test for a None response
        mock_get_compounds.return_value = []
        self.assertIsNone(name_to_smiles("unknown"))
if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTest(unittest.makeSuite(TestChemFunctions))
    runner = unittest.TextTestRunner()
    runner.run(suite)


#SMILESTOMOLARMASS
from ToolsKit.Project_ppchem import smiles_to_molar_mass

class TestMolecularFunctions(unittest.TestCase):
    @patch('Project_ppchem.Chem.MolFromSmiles')
    @patch('Project_ppchem.Descriptors.ExactMolWt')
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
    suite = unittest.TestSuite()
    suite.addTest(unittest.makeSuite(TestMolecularFunctions))
    runner = unittest.TextTestRunner()
    runner.run(suite)


#LINEARREG & MAKEGRAPH
from ToolsKit.Project_ppchem import linear_regression, make_graph

class TestPlottingFunctions(unittest.TestCase):
    @patch('Project_ppchem.pd.read_excel')
    @patch('Project_ppchem.plt')
    def test_linear_regression(self, mock_plt, mock_read_excel):
        # Setup
        mock_df = MagicMock()
        mock_df.values.tolist.return_value = [[1, 2], [3, 4]]
        mock_read_excel.return_value = mock_df

        # Test
        linear_regression("path.xlsx", "X", "Y", "Test")
        mock_plt.show.assert_called_once()

    @patch('Project_ppchem.pd.read_excel')
    @patch('Project_ppchem.plt')
    def test_make_graph(self, mock_plt, mock_read_excel):
        # Setup
        mock_df = MagicMock()
        mock_df.values.tolist.return_value = [[1, 2], [3, 4]]
        mock_read_excel.return_value = mock_df

        # Test
        make_graph("path.xlsx", "X", "Y", "Test")
        mock_plt.show.assert_called_once()
if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTest(unittest.makeSuite(TestPlottingFunctions))
    runner = unittest.TextTestRunner()
    runner.run(suite)
