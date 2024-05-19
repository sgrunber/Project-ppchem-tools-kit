selected_radio=None
from project_ppchem_tools_kit.clear_input import clear_input

def on_radio_select(value): 
    """Updates the selected_radio variable with the selected value and clears the content of the input field.

    Args:
        value (str): Value of the radio button that has been selected.
    """
    selected_radio.set(value)
    clear_input()