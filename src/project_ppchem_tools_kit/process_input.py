
from project_ppchem_tools_kit.name_to_smiles import name_to_smiles
from project_ppchem_tools_kit.linear_regression import linear_regression
from project_ppchem_tools_kit.smiles_to_molar_mass import smiles_to_molar_mass
from project_ppchem_tools_kit.make_graph import make_graph
from project_ppchem_tools_kit.display_molecule import display_molecule
from project_ppchem_tools_kit.error_calculation_interface import error_calculation_interface
from project_ppchem_tools_kit.display_result import display_result
entry_input = None

def process_input(event=None):
    """Processes the input based on the selected mode (molecule name, SMILES, file path for plotting). Validates the input and displays the appropriate results or errors.

    Args:
        event (Event, optional): The event that triggered this function, such as a button click or another event from a GUI component. Defaults to None.

    Returns:
        None: Does not return any values, directly affects the GUI by displaying messages or updating GUI components.
    """
    global selected_radio, x_label, y_label, title, result_text
    if selected_radio.get() == "1":
        smiles_molecule = name_to_smiles(entry_input.get().strip())
        if smiles_molecule is not None:
            result_text = f"The SMILES for the molecule '{entry_input.get().strip()}' is: {smiles_molecule}"
        else:
            result_text = f"The molecule '{entry_input.get().strip()}' was not found in the PubChem database."
    elif selected_radio.get() == "2":
        molar_mass = smiles_to_molar_mass(entry_input.get().strip())
        if molar_mass is not None:
            result_text = f"The molar mass of the molecule is: {molar_mass} g/mol"
        else:
            result_text = "Invalid SMILES or molecule not found."
    elif selected_radio.get() in ["3", "4"]:
        if selected_radio.get() == "3":
            make_graph(entry_input, x_label, y_label, title)
        else:
            linear_regression(entry_input, x_label, y_label, title)
        return  # Return to prevent displaying the result window
    elif selected_radio.get() == "5":  # Ajouter cette condition pour le bouton radio "5"
        display_molecule(entry_input)  # Appeler la fonction pour afficher la structure moléculaire
        return  # Return to prevent displaying the result window

    elif selected_radio.get() == "6":
        error_calculation_interface()
    else:
        result_text = "Please select an input type."

    display_result(result_text)