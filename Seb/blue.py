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


OUTPUT_PATH = Path(__file__).parent
ASSETS_PATH = OUTPUT_PATH / Path(r"/Users/gruenbergsebastien/Project-ppchem-tools-kit/tkinter/build/assets/frame0")


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
            result_text = f"The molar mass of the molecule is: {molar_mass}"
        else:
            result_text = "Invalid SMILES or molecule not found."
    elif selected_radio.get() == "3":
        file_path = entry_1.get().strip()  
        x_label = entry_3.get().strip()  
        y_label = entry_2.get().strip()  
        make_graph(file_path, x_label, y_label)
        #result_text = f"Graph created from file: {file_path}, X Label: {x_label}, Y Label: {y_label}"
    else:
        result_text = "Please select an input type."

    display_result(result_text)

def make_graph(filepath, x_label, y_label):
    try:
        df = pd.read_excel(filepath)
        data_list = df.values.tolist()
        x_values = [row[0] for row in data_list]  
        y_values = [row[1] for row in data_list]  
        plt.plot(x_values, y_values, 'k-')
        plt.xlabel(x_label)  
        plt.ylabel(y_label)  
        plt.title("Graph")
        plt.grid(True)
        plt.show()
    except FileNotFoundError:
        messagebox.showerror("Error", "The data file was not found.")


def display_result(result_text):
    result_window = tk.Toplevel(window)
    result_window.title("Result")

    result_textbox = tk.Text(result_window, wrap="word", font=("Times New Roman", 16), height=5, width=50)
    result_textbox.insert("1.0", result_text)
    result_textbox.config(state="disabled")
    result_textbox.grid(row=0, column=0, padx=50, pady=50)

    result_textbox.bind("<Control-a>", select_all)  
    result_textbox.bind("<Control-c>", copy_text)  

def browse_excel_file():
    filepath = filedialog.askopenfilename(title="Select Excel File", filetypes=(("Excel files", "*.xlsx"), ("All files", "*.*")))
    if filepath:
        try:
            entry_1.delete(0, tk.END)  # Effacer le contenu actuel de l'entry
            entry_1.insert(0, filepath)  # Insérer le chemin du fichier dans l'entry
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

input_type = tk.IntVar()
input_type.set(1)

def welcome_message():
    # Créer une nouvelle fenêtre pour afficher le message de bienvenue
    welcome_window = tk.Toplevel(window)
    welcome_window.title("Welcome Message")
    welcome_window.geometry("300x100")
    
    # Label contenant le message de bienvenue
    welcome_label = tk.Label(welcome_window, text="Welcome!", font=("Times New Roman", 40))
    welcome_label.pack(pady=20)

def bind_enter(event):
    process_input()

window.bind('<Return>', bind_enter)

window.geometry("1000x800")
window.configure(bg = "#75AAC8")


canvas = Canvas(
    window,
    bg = "#75AAC8",
    height = 800,
    width = 1000,
    bd = 0,
    highlightthickness = 0,
    relief = "ridge"
)

canvas.place(x = 0, y = 0)
canvas.create_text(
    12.0,
    177.0,
    anchor="nw",
    text="Choose Input Type :",
    fill="#FFFFFF",
    font=("Times New Roman", 27 * -1)
)



canvas.create_text(
    92.0,
    378.0,
    anchor="nw",
    text="Input :",
    fill="#FFFFFF",
    font=("Times New Roman", 30 * -1)
)

entry_image_1 = PhotoImage(
    file=relative_to_assets("entry_1.png"))
entry_bg_1 = canvas.create_image(
    572.5,
    401.5,
    image=entry_image_1
)
entry_1 = Entry(
    bd=3,
    bg="#8BD2EC",
    fg="white",
    highlightthickness=0,
    font=("Times New Roman", 25) 
)
entry_1.place(
    x=270.5,
    y=365.0,
    width=604.0,
    height=71.0
)

button_image_1 = PhotoImage(
    file=relative_to_assets("button_1.png"))
button_1 = Button(
    image=button_image_1,
    borderwidth=0,
    highlightthickness=0,
    command=browse_excel_file, 
    relief="flat"
)

button_1.place(
    x=186.0,
    y=528.0,
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
                                  font=("Times New Roman", 20), bg="#75AAC8")
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
    file=relative_to_assets("button_3.png"))
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
    383.0,
    500.0,
    anchor="nw",
    text="Y Axis Label :",
    fill="#FFFFFF",
    font=("Times New Roman", 20 * -1)
)

entry_image_2 = PhotoImage(
    file=relative_to_assets("entry_2.png"))
entry_bg_2 = canvas.create_image(
    657.0,
    520.5,
    image=entry_image_2
)
entry_2 = Entry(
    bd=1,
    bg="#8BD2EC",
    fg="white",
    highlightthickness=0,
    font=("Times New Roman", 20) 
)
entry_2.place(
    x=562.0,
    y=501.0,
    width=190.0,
    height=37.0
)

canvas.create_text(
    383.0,
    576.0,
    anchor="nw",
    text="X Axis Label :",
    fill="#FFFFFF",
    font=("Times New Roman", 20 * -1)
)

entry_image_3 = PhotoImage(
    file=relative_to_assets("entry_3.png"))
entry_bg_3 = canvas.create_image(
    657.0,
    596.5,
    image=entry_image_3
)
entry_3 = Entry(
    bd=1,
    bg="#8BD2EC",
    fg="white",
    highlightthickness=0,
    font=("Times New Roman", 20) 
)
entry_3.place(
    x=562.0,
    y=577.0,
    width=190.0,
    height=37.0
)



window.resizable(False, False)
window.mainloop()
