import threading
import time
import unittest
import tkinter as tk
from unittest.mock import patch

import sys
sys.path.insert(0, "./src")
from Chem_pack.error_calculation_interface import error_calculation_interface
from Chem_pack.error_calculation_interface import calculate_error_propagation, copy_latex_code, calculate_and_display, add_variable_entry, remove_variable_entry

class TestCalculateErrorPropagation(unittest.TestCase):

    def test_calculate_error_propagation(self):
        # Test data
        derivatives = [1, 2, 3]
        uncertainties = [0.1, 0.2, 0.3]
        
        # Expected result
        expected_error = 0.14907119849998599
        expected_average = 2.0
        
        # Calculate error propagation
        error, average = calculate_error_propagation(derivatives, uncertainties)
        
        # Check the result
        self.assertAlmostEqual(error, expected_error)
        self.assertAlmostEqual(average, expected_average)

class TestCopyLatexCode(unittest.TestCase):

    @patch('pyperclip.copy')
    def test_copy_latex_code(self, mock_copy):
        # Call the function
        copy_latex_code("Test LaTeX Code")

        # Assert that the copy function was called with the correct argument
        mock_copy.assert_called_once_with("Test LaTeX Code")

class TestCalculateAndDisplay(unittest.TestCase):

    # This is just a placeholder for now. You can expand it to include more comprehensive testing.
    def test_calculate_and_display(self):
        # Add your test cases here
        pass

class TestAddVariableEntry(unittest.TestCase):

    # This is just a placeholder for now. You can expand it to include more comprehensive testing.
    def test_add_variable_entry(self):
        # Add your test cases here
        pass

class TestRemoveVariableEntry(unittest.TestCase):

    # This is just a placeholder for now. You can expand it to include more comprehensive testing.
    def test_remove_variable_entry(self):
        # Add your test cases here
        pass

if __name__ == '__main__':
    unittest.main()


'''import unittest
from unittest.mock import patch, MagicMock
import pyperclip
import tkinter as tk

import sys
sys.path.insert(0, "./src")
from Chem_pack.error_calculation_interface import error_calculation_interface


class TestErrorCalculationInterface(unittest.TestCase):
    def setUp(self):
        # Patch the required tkinter functions and classes
        self.patcher_tk = patch('tkinter.Tk', return_value=MagicMock())
        self.mock_tk = self.patcher_tk.start()

        self.patcher_toplevel = patch('tkinter.Toplevel', return_value=MagicMock())
        self.mock_toplevel = self.patcher_toplevel.start()

        self.patcher_button = patch('tkinter.Button', return_value=MagicMock())
        self.mock_button = self.patcher_button.start()

        self.patcher_label = patch('tkinter.Label', return_value=MagicMock())
        self.mock_label = self.patcher_label.start()

        self.patcher_entry = patch('tkinter.Entry', return_value=MagicMock())
        self.mock_entry = self.patcher_entry.start()

        self.patcher_text = patch('tkinter.Text', return_value=MagicMock())
        self.mock_text = self.patcher_text.start()

        self.patcher_pyperclip = patch('pyperclip.copy', autospec=True)
        self.mock_pyperclip = self.patcher_pyperclip.start()

        # Initialize the lists as they are globals in the module
        global entries, entry_vars, entry_values, entry_uncertainties, remove_buttons, labels, error_calc_window, result_window
        entries = []
        entry_vars = []
        entry_values = []
        entry_uncertainties = []
        remove_buttons = []
        labels = []

        # Create the root window and set it as the global error_calc_window
        error_calc_window = self.mock_tk.return_value

        # Call the interface function to populate global variables
        error_calculation_interface()

    def tearDown(self):
        self.patcher_tk.stop()
        self.patcher_toplevel.stop()
        self.patcher_button.stop()
        self.patcher_label.stop()
        self.patcher_entry.stop()
        self.patcher_text.stop()
        self.patcher_pyperclip.stop()

    def test_calculate_and_display(self):
        # Define a local function to match the structure
        def calculate_and_display():
            data = []
            values = []
            uncertainties = []
            for i in range(len(entries)):
                var_name = entry_vars[i].get()
                var_value = float(entry_values[i].get())
                var_uncertainty = float(entry_uncertainties[i].get())
                data.append((var_name, var_value, var_uncertainty))
                values.append(var_value)
                uncertainties.append(var_uncertainty)

            # Calculate the propagation of uncertainty
            error, average = calculate_error_propagation(values, uncertainties)

            # Display the result with calculation details
            result_text = f"Propagation of Uncertainty (∆z): {error:.6f}\n"
            result_text += f"Average Value: {average:.6f}\n\n"
            result_text += "Calculation Details:\n"
            for var_name, var_value, var_uncertainty in data:
                result_text += f"{var_name}: Value: {var_value}, Uncertainty: {var_uncertainty}\n"

            latex_code = r"""
            Propagation of Uncertainty ($\Delta z$)
            The propagation of uncertainty is calculated as follows:
            \begin{equation}
            \Delta \bar{z} = \bar{z} \times \sqrt{\left(\left(\frac{\Delta x}{x}\right)^2 + \left(\frac{\Delta y}{y}\right)^2 + \ldots\right)}
            \end{equation}
            \[
            \Delta z = """ + f"{error:.6f}" + r"""
            \]

            The average value is calculated as follows:
            
            \begin{equation}
            \bar{z} = \frac{1}{N} \sum_{i=1}^{N} z_i
            \end{equation}
            \[
            \bar{z} = """ + f"{average:.6f}" + r"""
            \]

            The details of variables:
            \begin{itemize}
            """ + "\n".join([f"    \item {var_name}: Value: {var_value}, Uncertainty: {var_uncertainty}" for var_name, var_value, var_uncertainty in data]) + r"""
            \end{itemize}
            """

            # Create a new window to display the result
            global result_window
            result_window = tk.Toplevel(error_calc_window)
            result_window.title("Result")

            # Display the result in a Text widget
            result_textbox = tk.Text(result_window, wrap="word", font=("Times New Roman", 20), height=15, width=70)
            result_textbox.insert("1.0", result_text)
            result_textbox.grid(row=0, column=0, padx=5, pady=5)

            # Button to copy LaTeX code to clipboard
            copy_button = tk.Button(result_window, text="Copy LaTeX Code", command=lambda: copy_latex_code(latex_code))
            copy_button.grid(row=1, column=0, padx=5, pady=5)

        def calculate_error_propagation(derivatives, uncertainties):
            average_values = sum(derivatives) / len(derivatives)
            error = sum((unc / deriv) ** 2 for deriv, unc in zip(derivatives, uncertainties))
            standard_dev = (error ** 0.5) * average_values
            return standard_dev, average_values

        mock_entry_var = MagicMock()
        mock_entry_var.get.return_value = "Var1"
        mock_entry_value = MagicMock()
        mock_entry_value.get.return_value = "10.0"
        mock_entry_uncertainty = MagicMock()
        mock_entry_uncertainty.get.return_value = "0.5"
        
        entry_vars.append(mock_entry_var)
        entry_values.append(mock_entry_value)
        entry_uncertainties.append(mock_entry_uncertainty)
        entries.append((mock_entry_var, mock_entry_value, mock_entry_uncertainty))

        with patch('tkinter.Toplevel', return_value=MagicMock()) as mock_toplevel, \
             patch('tkinter.Text', return_value=MagicMock()) as mock_text, \
             patch('tkinter.Button', return_value=MagicMock()) as mock_button:
            
            calculate_and_display()
            mock_result_window = mock_toplevel.return_value
            mock_text.assert_called_once()
            mock_button.assert_called_once()
            
            self.assertTrue(mock_result_window.title.called)
            self.assertTrue(mock_text.called)
            self.assertTrue(mock_button.called)

    def test_remove_variable_entry(self):
        # Define a local function to match the structure
        def remove_variable_entry(row):
            entry = entries.pop(row - 1)
            for widget in entry:
                widget.destroy()
            for label in labels[(row - 1) * 3:(row - 1) * 3 + 3]:
                label.destroy()
            remove_buttons.pop(row - 1).destroy()

        mock_entry_var = MagicMock()
        mock_entry_value = MagicMock()
        mock_entry_uncertainty = MagicMock()
        
        entries.append((mock_entry_var, mock_entry_value, mock_entry_uncertainty))
        labels.extend([MagicMock(), MagicMock(), MagicMock()])
        remove_buttons.append(MagicMock())

        with patch.object(entries[0][0], 'destroy'), \
             patch.object(entries[0][1], 'destroy'), \
             patch.object(entries[0][2], 'destroy'), \
             patch.object(labels[0], 'destroy'), \
             patch.object(labels[1], 'destroy'), \
             patch.object(labels[2], 'destroy'), \
             patch.object(remove_buttons[0], 'destroy'):

            remove_variable_entry(1)

            self.assertEqual(len(entries), 0)
            self.assertEqual(len(labels), 0)
            self.assertEqual(len(remove_buttons), 0)

if __name__ == '__main__':
    unittest.main()

'''

'''
class TestErrorCalculationInterface(unittest.TestCase):
    def setUp(self):
        # Patch the required tkinter functions and classes
        self.patcher_tk = patch('tkinter.Tk', return_value=MagicMock())
        self.mock_tk = self.patcher_tk.start()

        self.patcher_toplevel = patch('tkinter.Toplevel', return_value=MagicMock())
        self.mock_toplevel = self.patcher_toplevel.start()

        self.patcher_button = patch('tkinter.Button', return_value=MagicMock())
        self.mock_button = self.patcher_button.start()

        self.patcher_label = patch('tkinter.Label', return_value=MagicMock())
        self.mock_label = self.patcher_label.start()

        self.patcher_entry = patch('tkinter.Entry', return_value=MagicMock())
        self.mock_entry = self.patcher_entry.start()

        self.patcher_text = patch('tkinter.Text', return_value=MagicMock())
        self.mock_text = self.patcher_text.start()

        self.patcher_pyperclip = patch('pyperclip.copy', autospec=True)
        self.mock_pyperclip = self.patcher_pyperclip.start()

        # Initialize the lists as they are globals in the module
        global entries, entry_vars, entry_values, entry_uncertainties, remove_buttons, labels, error_calc_window, result_window
        entries = []
        entry_vars = []
        entry_values = []
        entry_uncertainties = []
        remove_buttons = []
        labels = []

        # Create the root window and set it as the global error_calc_window
        error_calc_window = self.mock_tk.return_value

        # Call the interface function to populate global variables
        error_calculation_interface()

    def tearDown(self):
        self.patcher_tk.stop()
        self.patcher_toplevel.stop()
        self.patcher_button.stop()
        self.patcher_label.stop()
        self.patcher_entry.stop()
        self.patcher_text.stop()
        self.patcher_pyperclip.stop()

    def test_calculate_error_propagation(self):
        # Define a local function to match the structure
        def calculate_error_propagation(derivatives, uncertainties):
            average_values = sum(derivatives) / len(derivatives)
            error = sum((unc / deriv) ** 2 for deriv, unc in zip(derivatives, uncertainties))
            standard_dev = (error ** 0.5) * average_values
            return standard_dev, average_values

        derivatives = [1.0, 2.0, 3.0]
        uncertainties = [0.1, 0.2, 0.3]
        expected_standard_dev = (sum((u / d) ** 2 for d, u in zip(derivatives, uncertainties)) ** 0.5) * (sum(derivatives) / len(derivatives))
        expected_average_values = sum(derivatives) / len(derivatives)

        standard_dev, average_values = calculate_error_propagation(derivatives, uncertainties)

        self.assertAlmostEqual(standard_dev, expected_standard_dev)
        self.assertAlmostEqual(average_values, expected_average_values)

    def test_copy_latex_code(self):
        # Define a local function to match the structure
        def copy_latex_code(latex_code):
            pyperclip.copy(latex_code)
            result_window.clipboard_clear()
            result_window.clipboard_append(latex_code)
            result_window.update()

        latex_code = "Some LaTeX Code"
        result_window = MagicMock()

        with patch('tkinter.Toplevel', return_value=result_window):
            copy_latex_code(latex_code)

        self.mock_pyperclip.assert_called_once_with(latex_code)
        result_window.clipboard_clear.assert_called_once()
        result_window.clipboard_append.assert_called_once_with(latex_code)
        result_window.update.assert_called_once()

    def test_calculate_and_display(self):
        # Define a local function to match the structure
        def calculate_and_display():
            data = []
            values = []
            uncertainties = []
            for i in range(len(entries)):
                var_name = entry_vars[i].get()
                var_value = float(entry_values[i].get())
                var_uncertainty = float(entry_uncertainties[i].get())
                data.append((var_name, var_value, var_uncertainty))
                values.append(var_value)
                uncertainties.append(var_uncertainty)

            # Calculate the propagation of uncertainty
            error, average = calculate_error_propagation(values, uncertainties)

            # Display the result with calculation details
            result_text = f"Propagation of Uncertainty (∆z): {error:.6f}\n"
            result_text += f"Average Value: {average:.6f}\n\n"
            result_text += "Calculation Details:\n"
            for var_name, var_value, var_uncertainty in data:
                result_text += f"{var_name}: Value: {var_value}, Uncertainty: {var_uncertainty}\n"

            latex_code = r"""
            Propagation of Uncertainty ($\Delta z$)
            The propagation of uncertainty is calculated as follows:
            \begin{equation}
            \Delta \bar{z} = \bar{z} \times \sqrt{\left(\left(\frac{\Delta x}{x}\right)^2 + \left(\frac{\Delta y}{y}\right)^2 + \ldots\right)}
            \end{equation}
            \[
            \Delta z = """ + f"{error:.6f}" + r"""
            \]

            The average value is calculated as follows:
            
            \begin{equation}
            \bar{z} = \frac{1}{N} \sum_{i=1}^{N} z_i
            \end{equation}
            \[
            \bar{z} = """ + f"{average:.6f}" + r"""
            \]

            The details of variables:
            \begin{itemize}
            """ + "\n".join([f"    \item {var_name}: Value: {var_value}, Uncertainty: {var_uncertainty}" for var_name, var_value, var_uncertainty in data]) + r"""
            \end{itemize}
            """

            # Create a new window to display the result
            global result_window
            result_window = tk.Toplevel(error_calc_window)
            result_window.title("Result")

            # Display the result in a Text widget
            result_textbox = tk.Text(result_window, wrap="word", font=("Times New Roman", 20), height=15, width=70)
            result_textbox.insert("1.0", result_text)
            result_textbox.grid(row=0, column=0, padx=5, pady=5)

            # Button to copy LaTeX code to clipboard
            copy_button = tk.Button(result_window, text="Copy LaTeX Code", command=lambda: copy_latex_code(latex_code))
            copy_button.grid(row=1, column=0, padx=5, pady=5)

        def calculate_error_propagation(derivatives, uncertainties):
            average_values = sum(derivatives) / len(derivatives)
            error = sum((unc / deriv) ** 2 for deriv, unc in zip(derivatives, uncertainties))
            standard_dev = (error ** 0.5) * average_values
            return standard_dev, average_values

        mock_entry_var = MagicMock()
        mock_entry_var.get.return_value = "Var1"
        mock_entry_value = MagicMock()
        mock_entry_value.get.return_value = "10.0"
        mock_entry_uncertainty = MagicMock()
        mock_entry_uncertainty.get.return_value = "0.5"
        
        entry_vars.append(mock_entry_var)
        entry_values.append(mock_entry_value)
        entry_uncertainties.append(mock_entry_uncertainty)
        entries.append((mock_entry_var, mock_entry_value, mock_entry_uncertainty))

        with patch('tkinter.Toplevel', return_value=MagicMock()) as mock_toplevel:
            calculate_and_display()
            mock_result_window = mock_toplevel.return_value
            
            self.assertTrue(mock_result_window.title.called)
            self.assertTrue(mock_result_window.Text.called)
            self.assertTrue(mock_result_window.Button.called)

    def test_add_variable_entry(self):
        # Define a local function to match the structure
        def add_variable_entry():
            row = len(entries) + 1

            var_name_label = tk.Label(error_calc_window, text="Variable Name", bg="#1B262C", fg="white")
            var_name_label.grid(row=row, column=0, padx=5, pady=5)
            var_name_entry = tk.Entry(error_calc_window)
            var_name_entry.grid(row=row, column=1, padx=5, pady=5)
            entry_vars.append(var_name_entry)
            labels.append(var_name_label)

            var_value_label = tk.Label(error_calc_window, text="Value", bg="#1B262C", fg="white")
            var_value_label.grid(row=row, column=2, padx=5, pady=5)
            var_value_entry = tk.Entry(error_calc_window)
            var_value_entry.grid(row=row, column=3, padx=5, pady=5)
            entry_values.append(var_value_entry)
            labels.append(var_value_label)

            var_uncertainty_label = tk.Label(error_calc_window, text="Uncertainty", bg="#1B262C", fg="white")
            var_uncertainty_label.grid(row=row, column=4, padx=5, pady=5)
            var_uncertainty_entry = tk.Entry(error_calc_window)
            var_uncertainty_entry.grid(row=row, column=5, padx=5, pady=5)
            entry_uncertainties.append(var_uncertainty_entry)
            labels.append(var_uncertainty_label)

            entries.append((var_name_entry, var_value_entry, var_uncertainty_entry))

            remove_variable_button = tk.Button(error_calc_window, text="Remove", command=lambda r=row: remove_variable_entry(r),
                                                bg="#FF0000", fg="#1B262C")
            remove_variable_button.grid(row=row, column=6, padx=5, pady=5)
            remove_buttons.append(remove_variable_button)

        with patch('tkinter.Label', return_value=MagicMock()) as mock_label, \
             patch('tkinter.Entry', return_value=MagicMock()) as mock_entry, \
             patch('tkinter.Button', return_value=MagicMock()):

            add_variable_entry()

            self.assertTrue(mock_label.called)
            self.assertTrue(mock_entry.called)
            self.assertEqual(len(entries), 1)
            self.assertEqual(len(entry_vars), 1)
            self.assertEqual(len(entry_values), 1)
            self.assertEqual(len(entry_uncertainties), 1)

    def test_remove_variable_entry(self):
        # Define a local function to match the structure
        def remove_variable_entry(row):
            entry = entries.pop(row - 1)
            for widget in entry:
                widget.destroy()
            for label in labels[(row - 1) * 3:(row - 1) * 3 + 3]:
                label.destroy()
            remove_buttons.pop(row - 1).destroy()

        mock_entry_var = MagicMock()
        mock_entry_value = MagicMock()
        mock_entry_uncertainty = MagicMock()
        
        entries.append((mock_entry_var, mock_entry_value, mock_entry_uncertainty))
        labels.extend([MagicMock(), MagicMock(), MagicMock()])
        remove_buttons.append(MagicMock())

        remove_variable_entry(1)

        self.assertEqual(len(entries), 0)
        self.assertEqual(len(labels), 0)
        self.assertEqual(len(remove_buttons), 0)

if __name__ == '__main__':
    unittest.main()
'''