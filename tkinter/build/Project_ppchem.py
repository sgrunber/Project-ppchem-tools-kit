
import tkinter as tk
from sklearn.linear_model import LinearRegression
from sklearn.metrics import r2_score
from pathlib import Path
import numpy as np
from tkinter import ttk
from tkinter import messagebox
from tkinter import filedialog
import subprocess
import pubchempy as pcp
from rdkit import Chem
from rdkit.Chem import Descriptors
import matplotlib.pyplot as plt
import pandas as pd
from tkinter import Tk, Canvas, Entry, Text, Button, PhotoImage
from tkinter import IntVar, Radiobutton
from matplotlib_inline.backend_inline import set_matplotlib_formats
set_matplotlib_formats('png', 'pdf')

plt.rcParams['font.family'] = 'Times New Roman'

OUTPUT_PATH = Path(__file__).parent
ASSETS_PATH = OUTPUT_PATH / Path(r"/Users/meloenzinger/Documents/GitHub/Project-ppchem-tools-kit-bis/tkinter/build/assets/frame0")

def relative_to_assets(path: str) -> Path:
    """Constructs a path to a file located in the assets directory by combining the provided relative path with the ASSETS_PATH
    
    Args : 
        path (str) ; the relative path to be appended
    
    Returns : 
        the full path to the asset
    """
    return ASSETS_PATH / Path(path)


def name_to_smiles(molecule_name):
    """Fetches the SMILES representation of a molecule from its name from PubChem.

    Args:
        molecule_name (str): The name of the molecule.

    Returns:
        str: The canonical SMILES string if found, None otherwise.
    
    Raises:
        PubChemHTTPError: If any issue is encountered with the PubChem API request.
    """
    try:
        compound = pcp.get_compounds(molecule_name, 'name')
        if compound:
            return compound[0].canonical_smiles
        else:
            return None
    except pcp.PubChemHTTPError as e:
        print("Error occurred while fetching data from PubChem:", e)
        return None


def smiles_to_molar_mass(smiles):
    """Calculates the molar mass of a molecule from a given SMILES using RDKit.

    Args:
        smiles (str): The SMILES string of the molecule.

    Returns:
        float: The molar mass of the molecules (in grams per mole) if the SMILES string is valid, None otherwise.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is not None:
        return Descriptors.ExactMolWt(mol)
    else:
        return None



def linear_regression(file_path, x_label, y_label, title, grid=True, save_as=None, line_style='-', line_color='k'):
    """Performs a linear regression from an Excel file and plots the graph and R^2 value.

    Args:
        file_path (str): Path to the Excel file containing the data.
        x_label (str): x-axis label.
        y_label (_type_): y - axis label.
        title (str): Title of the graph.
        grid (bool): Whether to display grid lines on the plot. Defaults to True.
        save_as (str, optional): File path to save the graph. Defaults to None, in which case it is not saved.
        line_style (str, optional): Style of the line plot. Defaults to '-'.
        line_color (str, optional): Color of the plot line. Defaults to 'k' (black).
    
    Raises:
        FileNotFoundError: If the specified file_path does not exist
    """
    try:

        df = pd.read_excel(file_path)
        data_list = df.values.tolist()
        x_values = np.array([row[0] for row in data_list]).reshape(-1, 1) 
        y_values = [row[1] for row in data_list]

        model = LinearRegression()
        model.fit(x_values, y_values)


        y_pred = model.predict(x_values)

        r2 = r2_score(y_values, y_pred)

        plt.plot(x_values, y_pred, color='red', label='Linear Regression', linestyle=line_style, linewidth=1)
        plt.scatter(x_values, y_values, color='blue', label='Data Points')
        plt.xlabel(x_label, fontsize=15)
        plt.ylabel(y_label, fontsize=15)
        plt.title(title, fontsize=20)
        plt.rcParams['figure.dpi'] = 300
        plt.rcParams['savefig.dpi'] = 300
        plt.tight_layout()

        if grid:
            plt.grid(True)
        else:
            plt.grid(False)

        plt.legend()


        plt.text(0.6, 0.8, f'$R^2 = {r2:.2f}$', ha='center', va='center', transform=plt.gca().transAxes, fontsize=13, fontname='Times New Roman')

        plt.show()

        if save_as:
            plt.savefig(save_as)

    except FileNotFoundError:
        messagebox.showerror("Error", "The data file was not found.")


def make_graph(filepath, x_label, y_label, title, grid=True, save_as=None, line_style='-', line_color='k'):
    """Creates a scatter plot from data in a given Excel file.

    Args:
        filepath (str): Path to the Excel file containing the data.
        x_label (str): x-axis label.
        y_label (str): y- axis label.
        title (str): Plot title.
        grid (bool, optional): Whether to display grid lines on the plot. Defaults to True.
        save_as (str, optional): File path to save the plot. Defaults to None, in which case it is not saved.
        line_style (str, optional): Style of the line plot. Defaults to '-'.
        line_color (str, optional): Color of the plot line. Defaults to 'k' (black).
    """
    try:
       
        df = pd.read_excel(filepath)
        data_list = df.values.tolist()
        x_values = [row[0] for row in data_list]  
        y_values = [row[1] for row in data_list]  

        #plt.scatter(x_values, y_values, color='blue', label='Data Points')  (ici on pourait ajouté un truc qui permet d'affiché les point ou non)
        plt.plot(x_values, y_values, linestyle=line_style, color=line_color)
        plt.xlabel(x_label, fontsize = 15)  
        plt.ylabel(y_label, fontsize = 15) 
        plt.title(title, fontsize = 20)
        plt.rcParams['figure.dpi'] = 300
        plt.rcParams['savefig.dpi'] = 300
        plt.tight_layout()
        
    
        if grid:
            plt.grid(True)
        else:
            plt.grid(False)

        plt.show()

        if save_as:
            plt.savefig(save_as)

    except FileNotFoundError:
        messagebox.showerror("Error", "The data file was not found.")

def process_input(event=None):
    """Processes the input based on the selected mode (molecule name, SMILES, file path for plotting). Validates the input and displays the appropriate results or errors.

    Args:
        event (Event, optional): The event that triggered this function, such as a button click or another event from a GUI component. Defaults to None.

    Returns:
        None: Does not return any values, directly affects the GUI by displaying messages or updating GUI components.
    """
    input_text = entry_input.get().strip()
    if not input_text:
        messagebox.showerror("Error", "Please enter a molecule name, SMILES code, or file path.")
        return

    if selected_radio.get() == "1":
        smiles_molecule = name_to_smiles(input_text)
        if smiles_molecule is not None:
            result_text = f"The SMILES for the molecule '{input_text}' is: {smiles_molecule}"
        else:
            result_text = f"The molecule '{input_text}' was not found in the PubChem database."
    elif selected_radio.get() == "2":
        molar_mass = smiles_to_molar_mass(input_text)
        if molar_mass is not None:
            result_text = f"The molar mass of the molecule is: {molar_mass} g/mol"
        else:
            result_text = "Invalid SMILES or molecule not found."
    elif selected_radio.get() == "3":
        file_path = entry_input.get().strip()  
        x_label = entry_x_axis.get().strip()  
        y_label = entry_y_axis.get().strip()
        title = entry_graph_title.get().strip()  
        make_graph(file_path, x_label, y_label, title)
    elif selected_radio.get() == "4":  # Linear regression option
        file_path = entry_input.get().strip()  
        x_label = entry_x_axis.get().strip()  
        y_label = entry_y_axis.get().strip()
        title = entry_graph_title.get().strip()  
        linear_regression(file_path, x_label, y_label, title)
        return  # Return to prevent displaying the result window
    else:
        result_text = "Please select an input type."

    display_result(result_text)





def display_result(result_text):
    """Displays the result of an operation in a new window

    Args:
        result_text (str): The text to be displayed in the result window
    """
    result_window = tk.Toplevel(window)
    result_window.title("Result")


    result_textbox = tk.Text(result_window, wrap="word", font=("Times New Roman", 25), fg="#BBE1FA", bg = "#1B262C", height=10, width=30)
    result_textbox.insert("1.0", result_text)
    result_textbox.config(state="disabled")
    result_textbox.grid(row=0, column=0, padx=0, pady=0)




def browse_excel_file():
    """Opens a file dialog to select an Excel file, and updates the input field with the selected file path. Copies the file path to the clipboard.
    """
    filepath = filedialog.askopenfilename(title="Select Excel File", filetypes=(("Excel files", "*.xlsx"), ("All files", "*.*")))
    if filepath:
        try:
            entry_input.delete(0, tk.END)  
            entry_input.insert(0, filepath)  
            print("Selected Excel file:", filepath)
            print("File path copied to clipboard.")
        except Exception as e:
            print("Error:", e)

def select_all(event):
    """Selects all text in a widget, typically bound to a text-related widget event.

    Args:
        event: The event that triggered the function.

    Returns:
        None: Does not return any values.
    """
    event.widget.tag_add("sel", "1.0", "end")
    return "break"


def copy_text(event):
    """Copies selected text to the clipboard. Bound to an event.

    Args:
        event: The event that triggered the function.

    Returns:
        None: Does not return any values.
    """
    event.widget.event_generate("<<Copy>>")
    return "break"



window = tk.Tk()  #creates a Tkinter window instance

window.geometry(f"{window.winfo_reqwidth()}x{window.winfo_reqheight()}+{window.winfo_screenwidth()//2 - window.winfo_reqwidth()//2}+{window.winfo_screenheight()//2 - window.winfo_reqheight()//2}")
window.title("Project Tools-Kit")
input_type = tk.IntVar() #creates a Tkinter IntVar variable, which is used to track the value of the selected input type. In this case, it's initialized to 1
input_type.set(1)

def welcome_message():
    """Displays a welcome message in a top-level when triggered by a GUI event.
    """
    welcome_window = tk.Toplevel(window)
    welcome_window.title("Welcome Message")

    welcome_text = (
        "Welcome!\n\n"
        "Here is our project: https://github.com/sgrunber/Project-ppchem-tools-kit\n\n"
        "Enjoy ;)"
    )
    welcome_textbox = tk.Text(welcome_window, wrap="word", font=("Times New Roman", 25), fg="#BBE1FA", bg = "#1B262C", height=7, width=45)
    welcome_textbox.insert("1.0", welcome_text)
    welcome_textbox.config(state="disabled")
    welcome_textbox.grid(row=0, column=0, padx=0, pady=0)




window.geometry("1000x800")
window.configure(bg = "#1B262C")


canvas = Canvas( #creates an instance of the Canvas class in Tkinter
    window,
    bg = "#1B262C", #sets the background color of the Canvas
    height = 800, #sets the height of the Canvas in pixels
    width = 1000, #sets the width of the Canvas in pixels
    bd = 0, #sets the border size of the Canvas
    highlightthickness = 0, #sets the thickness of the Canvas highlight border
    relief = "ridge" #
)

canvas.place(x=0, y=0) #This places the Canvas at the coordinates (x, y)

canvas.create_text( #creates a text object on the Canvas
    12.0, # x coordinate
    194.0, # y coordinate
    anchor="nw",  #sets the anchor point of the text
    text="Choose Input Type :", #specifies the text to be displayed
    fill="#BBE1FA", #sets the fill color of the text
    font=("Times New Roman", 27, "bold") #specifies the font family, size, and weight of the text
)



canvas.create_text(  #principal input text label
    92.0,
    378.0,
    anchor="nw",    
    text="Input :",
    fill="#BBE1FA",
    font=("Times New Roman", 30, "bold")
)

entry_input = Entry( #principal input entry
    bd=1, 
    bg="#0F4C75",
    fg="#BBE1FA", #sets the text color of the Entry widget 
    highlightthickness=0,
    font=("Times New Roman", 25) 
)
entry_input.place( #principal input entry place
    x=272.5,
    y=368.0,
    width=605.0, #sets the width of the Entry widget in pixel
    height=74.0  #sets the height of the Entry widget in pixel
)



button_image_browse = PhotoImage( #browse button
    file=relative_to_assets("button_browse.png"))  
button_browse = Button(
    image=button_image_browse,
    borderwidth=0,
    highlightthickness=0,
    command=browse_excel_file, 
    relief="flat"
)

button_browse.place( #browse button place
    x=70.0,
    y=526.5,
    width=154.2904052734375,
    height=60.0
)

button_image_process = PhotoImage(  #Process button
    file=relative_to_assets("button_process.png"))
button_process = Button(
    image=button_image_process,
    borderwidth=0,
    highlightthickness=0,
    command=process_input,
    relief="flat"
)
button_process.place( #process button place
    x=409.0,
    y=680.0,
    width=182.0,
    height=73.0
)

def clear_input(): 
    """Clears the content of the input field in the GUI.
    """
    entry_input.delete(0, tk.END)

def on_radio_select(value): 
    """Updates the selected_radio variable with the selected value and clears the content of the input field.

    Args:
        value (str): Value of the radio button that has been selected.
    """
    selected_radio.set(value)
    clear_input()

def create_radio_button(x, y, text, value):
    """Creates a radio button with the specified text and value at the specified position (x, y)

    Args:
        x (int): The x-coordinate for placing the radio button.
        y (int): The y-coordinate for placing the radio button.
        text (str): The label text for the radio button.
        value (str): The value associated with the radio button.
    """
    radio_button = tk.Radiobutton(canvas, text=text, variable=selected_radio, value=value,
                                  command=lambda: on_radio_select(value),
                                  font=("Times New Roman", 20), bg="#1B262C")
    canvas.create_window(x, y, anchor="nw", window=radio_button)


radio_button_data = [
    (291, 194, "Molecule Name", "1"),
    (599, 194, "Excel Graph", "3"),
    (748, 194, "Linear Regression", "4"), 
    (774, 252, "Random", "random1"),
    (315, 252, "Random", "random2"),
    (469, 252, "Random", "random3"),
    (621, 252, "Random", "random4"),
    (469, 194, "SMILEs", "2")
]  
selected_radio = tk.StringVar(value="none")

for data in radio_button_data:
    create_radio_button(*data)


button_image_title = PhotoImage( #principale title button
    file=relative_to_assets("button_title.png"))
button_title = Button(
    image=button_image_title,
    borderwidth=0,
    highlightthickness=0,
    command=welcome_message, #command associated to this button
    relief="flat"
)
button_title.place( #principale title button place
    x=273.0,
    y=36.0,
    width=455.0,
    height=90.0
)


canvas.create_text(  #Y axis label text
    590.0,
    500.0,
    anchor="nw",
    text="Y Axis Label :",
    fill="#FFFFFF",
    font=("Times New Roman", 20 * -1)
)

entry_y_axis = Entry( #Y axis input
    bd=1,
    bg="#0F4C75",
    fg="#BBE1FA",
    highlightthickness=0,
    font=("Times New Roman", 20) 
)
entry_y_axis.place( # Y axis input place
    x=710.0,
    y=489.0,
    width=190.0,
    height=41.0
)


canvas.create_text( # Graph title label text
    300.0,
    570.0,
    anchor="nw",
    text="Title :",
    fill="#FFFFFF",
    font=("Times New Roman", 22 * -1)
)
entry_graph_title = Entry( #Graph title input
    bd=1,
    bg="#0F4C75",
    fg="#BBE1FA",
    highlightthickness=0,
    font=("Times New Roman", 20) 
)
entry_graph_title.place( #Graph title input place
    x=380.0,
    y=560.0,
    width=190.0,
    height=41.0
)


canvas.create_text(  #X axis label text
    260.0,
    500.0,
    anchor="nw",
    text="X Axis Label :",
    fill="#FFFFFF",
    font=("Times New Roman", 20 * -1)
)
entry_x_axis = Entry(  #X axis input
    bd=1,
    bg="#0F4C75",
    fg="#BBE1FA",
    highlightthickness=0,
    font=("Times New Roman", 20) 
)
entry_x_axis.place(  #X axis input place
    x=380.0,
    y=489.0,
    width=190.0,
    height=41.0
)


def bind_enter(event):
    """Binds the 'Enter' key to trigger the process_input function.

    Args:
        event: The event that triggered the function.
    """
    process_input()

window.bind('<Return>', bind_enter)


window.resizable(False, False)
window.mainloop()

