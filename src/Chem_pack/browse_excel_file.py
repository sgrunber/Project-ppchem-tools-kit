from tkinter import filedialog
import tkinter as tk

def browse_excel_file():
    """Opens a file dialog to select an Excel file, and updates the input field with the selected file path. Copies the file path to the clipboard.
    """
    global entry_input
    filepath = filedialog.askopenfilename(title="Select Excel File", filetypes=(("Excel files", "*.xlsx"), ("All files", "*.*")))
    if filepath:
        try:
            entry_input.delete(0, tk.END)  
            entry_input.insert(0, filepath)  
            print("Selected Excel file:", filepath)
            print("File path copied to clipboard.")
        except Exception as e:
            print("Error:", e)