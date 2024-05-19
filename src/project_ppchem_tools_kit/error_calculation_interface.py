# Chem_pack/error_calculation.py
import tkinter as tk
import pyperclip

def calculate_error_propagation(derivatives, uncertainties):
    average_values = sum(derivatives) / len(derivatives)
    error = sum((unc / deriv) ** 2 for deriv, unc in zip(derivatives, uncertainties))
    standard_dev = (error ** 0.5) * average_values
    return standard_dev, average_values

def error_calculation_interface():
    def copy_latex_code(latex_code):
        pyperclip.copy(latex_code)
        result_window.clipboard_clear()
        result_window.clipboard_append(latex_code)
        result_window.update()

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

        \subsection*{Propagation of Uncertainty ($\Delta\bar{z}$)}

        The propagation of uncertainty is calculated as follows:

        \begin{equation}
        \Delta \bar{z} = \bar{z} \times \sqrt{\left(\left(\frac{\Delta x}{x}\right)^2 + \left(\frac{\Delta y}{y}\right)^2 + \ldots\right)}
        \end{equation}

        The calculated uncertainty is:

        \[
        \Delta \bar{z} = """ + f"{error:.6f}" + r"""
        \]

        The average value is calculated as follows:

        \begin{equation}
        \bar{z} = \frac{1}{N} \sum_{i=1}^{N} z_i
        \end{equation}

        The calculated average value is:

        \[
        \bar{z} = """ + f"{average:.6f}" + r"""
        \]

        \subsection*{Details of Variables}

        \begin{itemize}
        """ + "\n".join([f"    \\item {var_name}: Value: {var_value}, Uncertainty: {var_uncertainty}" for var_name, var_value, var_uncertainty in data]) + r"""
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

    def add_variable_entry():
        # Ajouter des champs d'entrée pour une nouvelle variable
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

    def remove_variable_entry(row):
        entry = entries.pop(row - 1)
        for widget in entry:
            widget.destroy()
        for label in labels[(row - 1) * 3:(row - 1) * 3 + 3]:
            label.destroy()
        remove_buttons.pop(row - 1).destroy()

    error_calc_window = tk.Tk()
    error_calc_window.title("Calcul d'erreur")
    error_calc_window.configure(bg="#1B262C")

    add_variable_button = tk.Button(error_calc_window, text="Ajouter une variable", command=add_variable_entry, height=2,
                                    bg="#0000FF", fg="#1B262C", font=("Times New Roman", 15))
    add_variable_button.grid(row=0, column=2, columnspan=3, padx=5, pady=5, sticky="")

    # Pour le bouton "Calculer la propagation d'erreur"
    calculate_error_button = tk.Button(error_calc_window, text="Calculer la propagation d'erreur",
                                        command=calculate_and_display, height=2, bg="#008000", fg="#1B262C",
                                        font=("Times New Roman", 15))
    calculate_error_button.grid(row=20, column=2, columnspan=3, padx=5, pady=5, sticky="")

    # Listes pour stocker les widgets d'entrée
    entries = []
    entry_vars = []
    entry_values = []
    entry_uncertainties = []
    remove_buttons = []
    labels = []

    error_calc_window.mainloop()
