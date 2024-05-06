
import tkinter as tk
from pathlib import Path
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
from IPython.display import set_matplotlib_formats
set_matplotlib_formats('svg', 'pdf')

plt.rcParams['font.family'] = 'Times New Roman'

OUTPUT_PATH = Path(__file__).parent
ASSETS_PATH = OUTPUT_PATH / Path(r"/Users/meloenzinger/Desktop/EPFL/BA4/prog/Project-ppchem-tools-kit/tkinter/build/assets/frame0")

def relative_to_assets(path: str) -> Path:
    return ASSETS_PATH / Path(path)


def select_mode(mode):
    print("Selected mode:", mode)

def name_to_smiles(molecule_name):
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
    mol = Chem.MolFromSmiles(smiles)
    if mol is not None:
        return Descriptors.ExactMolWt(mol)
    else:
        return None


def process_input(event=None):
    input_text = entry_1.get().strip()
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
        file_path = entry_1.get().strip()  
        x_label = entry_3.get().strip()  
        y_label = entry_2.get().strip()
        title = entry_4.get().strip()  
        make_graph(file_path, x_label, y_label, title)
        
    else:
        result_text = "Please select an input type."

    display_result(result_text)



def make_graph(filepath, x_label, y_label, title, grid=True, save_as=None, line_style='-', line_color='k'):
    try:
       
        df = pd.read_excel(filepath)
        data_list = df.values.tolist()
        x_values = [row[0] for row in data_list]  
        y_values = [row[1] for row in data_list]  

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

def display_result(result_text):
    result_window = tk.Toplevel(window)
    result_window.title("Result")


    result_textbox = tk.Text(result_window, wrap="word", font=("Times New Roman", 25), fg="#BBE1FA", bg = "#1B262C", height=10, width=30)
    result_textbox.insert("1.0", result_text)
    result_textbox.config(state="disabled")
    result_textbox.grid(row=0, column=0, padx=0, pady=0)

    result_textbox.bind("<Control-a>", select_all)  
    result_textbox.bind("<Control-c>", copy_text)  

def browse_excel_file():
    filepath = filedialog.askopenfilename(title="Select Excel File", filetypes=(("Excel files", "*.xlsx"), ("All files", "*.*")))
    if filepath:
        try:
            entry_1.delete(0, tk.END)  
            entry_1.insert(0, filepath)  
            print("Selected Excel file:", filepath)
            print("File path copied to clipboard.")
        except Exception as e:
            print("Error:", e)

def select_all(event):
    event.widget.tag_add("sel", "1.0", "end")
    return "break"


def copy_text(event):
    event.widget.event_generate("<<Copy>>")
    return "break"



window = tk.Tk()

window.geometry(f"{window.winfo_reqwidth()}x{window.winfo_reqheight()}+{window.winfo_screenwidth()//2 - window.winfo_reqwidth()//2}+{window.winfo_screenheight()//2 - window.winfo_reqheight()//2}")

input_type = tk.IntVar()
input_type.set(1)

def welcome_message():
    welcome_window = tk.Toplevel(window)
    welcome_window.title("Welcome Message")
    welcome_window.geometry("600x200")
    
    
    welcome_label = tk.Label(welcome_window, text="Welcome! \n  \n Here is our project : https://github.com/sgrunber/Project-ppchem-tools-kit \n \n Enjoy ;)", font=("Times New Roman", 20)) 
    welcome_label.pack(pady=20)



def bind_enter(event):
    process_input()

window.bind('<Return>', bind_enter)

window.geometry("1000x800")
window.configure(bg = "#228B22")

canvas = Canvas(
    window,
    bg = "#1B262C",
    height = 800,
    width = 1000,
    bd = 0,
    highlightthickness = 0,
    relief = "ridge"
)

canvas.place(x=0, y=0)
canvas.create_text(
    12.0,
    194.0,
    anchor="nw",
    text="Choose Input Type :",
    fill="#BBE1FA",
    font=("Times New Roman", 27, "bold")  
)



canvas.create_text(
    92.0,
    378.0,
    anchor="nw",
    text="Input :",
    fill="#BBE1FA",
    font=("Times New Roman", 30, "bold")
)
'''
entry_image_1 = PhotoImage(
    file=relative_to_assets("entry_1.png"))
entry_bg_1 = canvas.create_image(
    574.0,
    402.5,
    image=entry_image_1
)
'''
entry_1 = Entry(
    bd=1,
    bg="#0F4C75",
    fg="#BBE1FA",
    highlightthickness=0,
    font=("Times New Roman", 25) 
)
entry_1.place(
    x=272.5,
    y=368.0,
    width=605.0,
    height=74.0
)

button_image_1 = PhotoImage(
    file=relative_to_assets("button_3.png"))
button_1 = Button(
    image=button_image_1,
    borderwidth=0,
    highlightthickness=0,
    command=browse_excel_file, 
    relief="flat"
)

button_1.place(
    x=70.0,
    y=526.5,
    width=154.2904052734375,
    height=60.0
)

button_image_2 = PhotoImage(
    file=relative_to_assets("button_2.png"))
button_2 = Button(
    image=button_image_2,
    borderwidth=0,
    highlightthickness=0,
    command=process_input,
    relief="flat"
)
button_2.place(
    x=409.0,
    y=680.0,
    width=182.0,
    height=73.0
)

def clear_input():
    entry_1.delete(0, tk.END)

def on_radio_select(value):
    selected_radio.set(value)
    clear_input()



# Fonction pour créer les Radiobuttons
def create_radio_button(x, y, text, value):
    radio_button = tk.Radiobutton(canvas, text=text, variable=selected_radio, value=value,
                                  command=lambda: on_radio_select(value),
                                  font=("Times New Roman", 20), bg="#1B262C")
    canvas.create_window(x, y, anchor="nw", window=radio_button)

# Coordonnées et textes pour les Radiobuttons
radio_button_data = [
    (291, 193, "Molecule Name", "1"),
    (599, 192, "Excel Graph", "3"),
    (748, 192, "Linear Regression", "linear"),
    (774, 252, "Random", "random1"),
    (315, 252, "Random", "random2"),
    (469, 252, "Random", "random3"),
    (621, 252, "Random", "random4"),
    (469, 192, "SMILEs", "2")
]  
# Variable pour suivre le Radiobutton sélectionné
selected_radio = tk.StringVar(value="none")

# Création des Radiobuttons à partir des données
for data in radio_button_data:
    create_radio_button(*data)

button_image_3 = PhotoImage(
    file=relative_to_assets("button_1.png"))
button_3 = Button(
    image=button_image_3,
    borderwidth=0,
    highlightthickness=0,
    command=welcome_message,
    relief="flat"
)
button_3.place(
    x=273.0,
    y=36.0,
    width=455.0,
    height=90.0
)

canvas.create_text(
    590.0,
    500.0,
    anchor="nw",
    text="Y Axis Label :",
    fill="#FFFFFF",
    font=("Times New Roman", 20 * -1)
)
'''
entry_image_2 = PhotoImage(
    file=relative_to_assets("entry_2.png"))
entry_bg_2 = canvas.create_image(
    657.0,
    520.5,
    image=entry_image_2
)
'''
entry_2 = Entry(
    bd=1,
    bg="#0F4C75",
    fg="#BBE1FA",
    highlightthickness=0,
    font=("Times New Roman", 20) 
)
entry_2.place(
    x=710.0,
    y=489.0,
    width=190.0,
    height=41.0
)

canvas.create_text(
    260.0,
    500.0,
    anchor="nw",
    text="X Axis Label :",
    fill="#FFFFFF",
    font=("Times New Roman", 20 * -1)
)
canvas.create_text(
    300.0,
    570.0,
    anchor="nw",
    text="Title :",
    fill="#FFFFFF",
    font=("Times New Roman", 22 * -1)
)
entry_4 = Entry(
    bd=1,
    bg="#0F4C75",
    fg="#BBE1FA",
    highlightthickness=0,
    font=("Times New Roman", 20) 
)
entry_4.place(
    x=380.0,
    y=560.0,
    width=190.0,
    height=41.0
)
'''
entry_image_3 = PhotoImage(
    file=relative_to_assets("entry_3.png"))
entry_bg_3 = canvas.create_image(
    657.0,
    596.5,
    image=entry_image_3
)
'''
entry_3 = Entry(
    bd=1,
    bg="#0F4C75",
    fg="#BBE1FA",
    highlightthickness=0,
    font=("Times New Roman", 20) 
)
entry_3.place(
    x=380.0,
    y=489.0,
    width=190.0,
    height=41.0
)





window.resizable(False, False)
window.mainloop()


