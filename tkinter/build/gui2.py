import tkinter as tk
from tkinter import ttk
from tkinter import messagebox
from tkinter import filedialog
import subprocess
import pubchempy as pcp
from rdkit import Chem
from rdkit.Chem import Descriptors
import matplotlib.pyplot as plt
import pandas as pd
from pathlib import Path
from tkinter import Tk, Canvas, Entry, Text, Button, PhotoImage


OUTPUT_PATH = Path(__file__).parent
ASSETS_PATH = OUTPUT_PATH / Path(r"/Users/gruenbergsebastien/Project-ppchem-tools-kit/tkinter/build/assets/frame0")


def relative_to_assets(path: str) -> Path:
    return ASSETS_PATH / Path(path)



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
    input_text = input_entry.get().strip()
    if not input_text:
        messagebox.showerror("Error", "Please enter a molecule name, SMILES code, or file path.")
        return

    if input_type.get() == 1:
        smiles_molecule = name_to_smiles(input_text)
        if smiles_molecule is not None:
            result_text = f"The SMILES for the molecule '{input_text}' is: {smiles_molecule}"
        else:
            result_text = f"The molecule '{input_text}' was not found in the PubChem database."
    elif input_type.get() == 2:
        molar_mass = smiles_to_molar_mass(input_text)
        if molar_mass is not None:
            result_text = f"The molar mass of the molecule is: {molar_mass}"
        else:
            result_text = "Invalid SMILES or molecule not found."
    elif input_type.get() == 3:
        file_path = entry_1.get().strip()  # Récupérer le chemin d'accès au fichier Excel
        x_label = x_label_entry.get().strip()  # Récupérer le nom de l'axe des abscisses
        y_label = y_label_entry.get().strip()  # Récupérer le nom de l'axe des ordonnées
        make_graph(file_path, x_label, y_label)  # Envoyer le chemin d'accès au fichier Excel à la fonction make_graph
        return
    else:
        result_text = "Please select an input type."

    display_result(result_text)


def make_graph(filepath, x_label, y_label):
    try:
        df = pd.read_excel(filepath)
        data_list = df.values.tolist()
        x_values = [row[0] for row in data_list]  # Première colonne
        y_values = [row[1] for row in data_list]  # Deuxième colonne
        plt.plot(x_values, y_values, 'k-')
        plt.xlabel(x_label)  # Utiliser le nom de l'axe des abscisses fourni
        plt.ylabel(y_label)  # Utiliser le nom de l'axe des ordonnées fourni
        plt.title("Graph")
        plt.grid(True)
        plt.show()
    except FileNotFoundError:
        messagebox.showerror("Error", "The data file was not found.")


def display_result(result_text):
    result_window = tk.Toplevel(root)
    result_window.title("Result")

    result_textbox = tk.Text(result_window, wrap="word", font=("Times New Roman", 16), height=5, width=50)
    result_textbox.insert("1.0", result_text)
    result_textbox.config(state="disabled")
    result_textbox.grid(row=0, column=0, padx=50, pady=50)

    result_textbox.bind("<Control-a>", select_all)  # Ctrl+A pour tout sélectionner
    result_textbox.bind("<Control-c>", copy_text)  # Ctrl+C pour copier


def browse_excel_file():
    filepath = filedialog.askopenfilename(title="Select Excel File", filetypes=(("Excel files", "*.xlsx"), ("All files", "*.*")))
    if filepath:
        try:
            # Exécuter la commande pour copier le chemin dans le presse-papiers
            subprocess.run(['pbcopy'], universal_newlines=True, input=filepath.strip())
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


def select_mode():
    mode = input_type.get()  

    if mode == 1:
        process = name_to_smiles
    elif mode == 2:
        process = smiles_to_molar_mass
    elif mode == 3:
        process = make_graph
    else:
        messagebox.showerror("Error", "Please select an input type.")
        return

    input_text = input_entry.get().strip()

    if not input_text:
        messagebox.showerror("Error", "Please enter a molecule name, SMILES code, or file path.")
        return

    process_input(process, input_text)

def process_input(process, input_text):
    result_text = process(input_text)
    display_result(result_text)



window = Tk()

window.geometry("1000x800")
window.configure(bg = "#75AAC8")

input_type = tk.IntVar()
input_type.set(1)
input_frame = tk.Frame(window)


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
    192.0,
    anchor="nw",
    text="Choose Input Type :",
    fill="#FFFFFF",
    font=("Times New Roman", 27 * -1)
)

button_image_1 = PhotoImage(
    file=relative_to_assets("button_1.png"))
button_1 = Button(
    image=button_image_1,
    borderwidth=0,
    highlightthickness=0,
    command=select_mode,
    relief="flat"
)
button_1.place(
    x=469.0,
    y=192.0,
    width=101.0,
    height=37.0
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
    bd=1,
    bg="#A9D7E6",
    fg="white", 
    highlightthickness=0,
    font=("Times New Roman", 25) 
)


entry_1.place(
    x=270.5 + 5,  
    y=365.0 + 5, 
    width=604.0 - 10,  
    height=71.0 - 10  
)
entry_1.insert(0, "Your text here")


button_image_2 = PhotoImage(
    file=relative_to_assets("button_2.png"))
button_2 = Button(
    image=button_image_2,
    borderwidth=0,
    highlightthickness=0,
    command=select_mode,
    relief="flat"
    )

button_2.place(
    x=413.0,
    y=513.0,
    width=182.0,
    height=73.0
)

button_image_3 = PhotoImage(
    file=relative_to_assets("button_3.png"))
button_3 = Button(
    image=button_image_3,
    borderwidth=0,
    highlightthickness=0,
    command=select_mode,
    relief="flat"
)


button_3.place(
    x=413.0,
    y=652.0,
    width=182.0,
    height=73.0
)

button_image_4 = PhotoImage(
    file=relative_to_assets("button_4.png"))
button_4 = Button(
    image=button_image_4,
    borderwidth=0,
    highlightthickness=0,
    command=select_mode,
    relief="flat"
)
button_4.place(
    x=289.0,
    y=192.0,
    width=153.0,
    height=37.0
)

button_image_5 = PhotoImage(
    file=relative_to_assets("button_5.png"))
button_5 = Button(
    image=button_image_5,
    borderwidth=0,
    highlightthickness=0,
    command=select_mode,
    relief="flat"
)
button_5.place(
    x=599.0,
    y=192.0,
    width=125.0,
    height=37.0
)

button_image_6 = PhotoImage(
    file=relative_to_assets("button_6.png"))
button_6 = Button(
    image=button_image_6,
    borderwidth=0,
    highlightthickness=0,
    command=select_mode,
    relief="flat"
)
button_6.place(
    x=748.0,
    y=192.0,
    width=175.0,
    height=37.0
)

button_image_7 = PhotoImage(
    file=relative_to_assets("button_7.png"))
button_7 = Button(
    image=button_image_7,
    borderwidth=0,
    highlightthickness=0,
    command=select_mode,
    relief="flat"
)
button_7.place(
    x=774.0,
    y=252.0,
    width=101.0,
    height=37.0
)

button_image_8 = PhotoImage(
    file=relative_to_assets("button_8.png"))
button_8 = Button(
    image=button_image_8,
    borderwidth=0,
    highlightthickness=0,
    command=select_mode,
    relief="flat"
)
button_8.place(
    x=315.0,
    y=252.0,
    width=101.0,
    height=37.0
)

button_image_9 = PhotoImage(
    file=relative_to_assets("button_9.png"))
button_9 = Button(
    image=button_image_9,
    borderwidth=0,
    highlightthickness=0,
    command=select_mode,
    relief="flat"
)
button_9.place(
    x=469.0,
    y=252.0,
    width=101.0,
    height=37.0
)

button_image_10 = PhotoImage(
    file=relative_to_assets("button_10.png"))
button_10 = Button(
    image=button_image_10,
    borderwidth=0,
    highlightthickness=0,
    command=browse_excel_file,
    relief="flat"
)
button_10.place(
    x=621.0,
    y=252.0,
    width=101.0,
    height=37.0
)

button_image_11 = PhotoImage(
    file=relative_to_assets("button_11.png"))
button_11 = Button(
    image=button_image_11,
    borderwidth=0,
    highlightthickness=0,
    command=process_input,
    relief="flat"
)
button_11.place(
    x=284.0,
    y=40.0,
    width=455.0,
    height=90.0
)


'''
molecule_name_radio.grid(row=0, column=1, padx=5, sticky="w")

smiles_radio = tk.Radiobutton(input_frame, text="SMILES", variable=input_type, value=2,
                               font=("Times New Roman", 20))
smiles_radio.grid(row=0, column=2, padx=5, sticky="w")

excel_radio = tk.Radiobutton(input_frame, text="Excel Data", variable=input_type, value=3,
                              font=("Times New Roman", 20))
excel_radio.grid(row=0, column=3, padx=5, sticky="w")

browse_button = ttk.Button(text="Browse Excel File", command=browse_excel_file, style="browse.TButton")
'''


window.resizable(False, False)
window.mainloop()
