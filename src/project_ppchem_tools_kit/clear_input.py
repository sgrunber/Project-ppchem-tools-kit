entry_input = None
import tkinter as tk

def clear_input(): 
    """Clears the content of the input field in the GUI.
    """
    entry_input.delete(0, tk.END)
    