import unittest
from unittest.mock import patch, MagicMock
import tkinter as tk

#importing functions from the project file
from Project_ppchem import relative_to_assets, ASSETS_PATH, name_to_smiles, smiles_to_molar_mass, linear_regression, make_graph, process_input, display_result, browse_excel_file, select_all, copy_text, welcome_message, clear_input, on_radio_select, create_radio_button, bind_enter
import Project_ppchem


#1.RELATIVETOASSETS
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


#2.NAMETOSMILES
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


#3.SMILESTOMOLARMASS
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


#4,5 LINEARREG & MAKEGRAPH
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
        with patch('Project_ppchem.entry_input', self.mock_entry_instance):
            process_input()
            self.mock_messagebox.showerror.assert_called_once_with("Error", "Please enter a molecule name, SMILES code, or file path.")
if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTest(unittest.makeSuite(TestProcessInput))
    runner = unittest.TextTestRunner()
    runner.run(suite)


#7. DISPLAYRESULT !!!REVIEW
class TestDisplayResult(unittest.TestCase):
    def setUp(self):
        # Patching the entire Entry class might be necessary if it's used like this throughout your app
        self.patcher_entry = patch('tkinter.Entry')
        self.mock_entry = self.patcher_entry.start()
        self.mock_entry_instance = MagicMock()
        self.mock_entry.return_value = self.mock_entry_instance
        self.mock_entry_instance.get.return_value = ''

        # Mock the messagebox to prevent actual GUI dialog boxes during the test
        self.patcher_msgbox = patch('tkinter.messagebox', autospec=True)
        self.mock_messagebox = self.patcher_msgbox.start()

        # Assuming 'selected_radio' is managed similarly, often as a part of a class instance
        self.patcher_radio = patch('tkinter.IntVar', autospec=True)
        self.mock_radio = self.patcher_radio.start()
        self.mock_radio_instance = MagicMock()
        self.mock_radio.return_value = self.mock_radio_instance
        self.mock_radio_instance.get.return_value = '1'  # Example value

    def tearDown(self):
        self.patcher_entry.stop()
        self.patcher_msgbox.stop()
        self.patcher_radio.stop()

    def test_empty_input(self):
        """Test process_input with empty input."""
        process_input()
        self.mock_messagebox.showerror.assert_called_once_with("Error", "Please enter a molecule name, SMILES code, or file path.")
if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTest(unittest.makeSuite(TestDisplayResult))
    runner = unittest.TextTestRunner()
    runner.run(suite)



#8. BROWSEEXCELFILE !!!REVIEW
import tkinter as tk
class TestBrowseExcelFile(unittest.TestCase):
    def test_browse_excel_file(self):
        """Test that browse_excel_file sets the entry with the file path."""
        with patch('tkinter.filedialog.askopenfilename', return_value='/path/to/file.xlsx') as mock_filedialog, \
             patch('tkinter.Entry', MagicMock()) as mock_entry:
            browse_excel_file()
            mock_entry.return_value.insert.assert_called_with(0, '/path/to/file.xlsx')
if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTest(unittest.makeSuite(TestBrowseExcelFile))
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


#11. WELCOMEMESSAGE !!!REVIEW ALL BELOW
import tkinter as tk
class TestWelcomeMessage(unittest.TestCase):
    @patch('tkinter.Toplevel')
    @patch('tkinter.Text')
    def test_welcome_message(self, mock_text, mock_toplevel):
        # Create a mock Tk root window
        mock_root = MagicMock(spec=tk.Tk)
        
        # Call the function that should use this mock
        welcome_message(mock_root)

        # Verify that Toplevel was called with the mock root window
        mock_toplevel.assert_called_once_with(mock_root) 
if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTest(unittest.makeSuite(TestWelcomeMessage))
    runner = unittest.TextTestRunner()
    runner.run(suite)


#12. CLEARINPUT    
import tkinter as tk
class TestClearInput(unittest.TestCase):
    @patch('tkinter.Entry', autospec=True)
    def test_clear_input(self, mock_entry):
        mock_entry_instance = MagicMock()
        mock_entry.return_value = mock_entry_instance
        
        clear_input()
        
        mock_entry_instance.delete.assert_called_once_with(0, tk.END)
if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTest(unittest.makeSuite(TestClearInput))
    runner = unittest.TextTestRunner()
    runner.run(suite)

#13. ONRADIOSELECT
import tkinter as tk
class TestOnRadioSelect(unittest.TestCase):
    @patch('tkinter.IntVar', autospec=True)
    @patch('tkinter.Entry', autospec=True)
    def test_on_radio_select(self, mock_entry, mock_intvar):
        mock_entry_instance = MagicMock()
        mock_entry.return_value = mock_entry_instance
        
        mock_intvar_instance = MagicMock()
        mock_intvar.return_value = mock_intvar_instance
        
        value = 'test_value'
        on_radio_select(value)
        
        mock_intvar_instance.set.assert_called_once_with(value)
        mock_entry_instance.delete.assert_called_once_with(0, tk.END)
if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTest(unittest.makeSuite(TestOnRadioSelect))
    runner = unittest.TextTestRunner()
    runner.run(suite)

#14. RADIOBUTTON
import tkinter as tk
class TestCreateRadioButton(unittest.TestCase):
    @patch('tkinter.Canvas', autospec=True)
    @patch('tkinter.Radiobutton', autospec=True)
    def test_create_radio_button(self, mock_radiobutton, mock_canvas):
        mock_canvas_instance = MagicMock()
        mock_canvas.return_value = mock_canvas_instance
        
        x, y = 100, 200
        text = 'Test Radio'
        value = 'radio_value'
        
        create_radio_button(x, y, text, value)
        
        mock_radiobutton.assert_called_once_with(mock_canvas_instance, text=text, variable=selected_radio, value=value,
                                                 command=lambda: on_radio_select(value),
                                                 font=("Times New Roman", 20), bg="#1B262C")
        mock_canvas_instance.create_window.assert_called_once_with(x, y, anchor="nw", window=mock_radiobutton.return_value)
if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTest(unittest.makeSuite(TestCreateRadioButton))
    runner = unittest.TextTestRunner()
    runner.run(suite)


#15. BINDENTER
import tkinter as tk
class TestBindEnter(unittest.TestCase):
    @patch('module_where_defined.process_input')
    @patch('tkinter.Tk')
    def test_bind_enter(self, mock_tk, mock_process_input):
        # Create a mock window instance from the patched tkinter Tk
        mock_window = MagicMock()
        mock_tk.return_value = mock_window
        
        # Setting up the bind
        window = mock_tk()
        window.bind('<Return>', bind_enter)
        
        # Simulate the Enter key event
        event = MagicMock()
        bind_enter(event)
        
        # Verify process_input was called
        mock_process_input.assert_called_once()
        
        # Check that the bind was set correctly
        mock_window.bind.assert_called_once_with('<Return>', bind_enter)
if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTest(unittest.makeSuite(TestBindEnter))
    runner = unittest.TextTestRunner()
    runner.run(suite)
