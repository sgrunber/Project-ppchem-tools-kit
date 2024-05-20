from rdkit import Chem
from rdkit.Chem import Draw
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
import tkinter as tk



def display_molecule(entry_input):
    """Displays the 2D structure of a molecule from a given SMILES string.

    Args:
        entry_input (tk.Entry): The entry widget containing the SMILES string.
    """
    smiles = entry_input.get().strip()  
    mol = Chem.MolFromSmiles(smiles)  

    if mol is not None:
        def show_molecule_2d():
            img = Draw.MolToImage(mol)

            mol_window = tk.Toplevel()
            mol_window.title("Molecular Structure")

            pimg = FigureCanvasTkAgg(plt.Figure(figsize=(4, 3)), master=mol_window)
            ax = pimg.figure.add_subplot(111)
            ax.imshow(img, interpolation='bilinear')
            ax.axis('off')

            canvas = pimg.get_tk_widget()
            canvas.pack()

            toolbar = NavigationToolbar2Tk(pimg, mol_window)  
            toolbar.update()
            toolbar.pack()

            mol_window.mainloop()

        show_molecule_2d()
    else:
        print("Error : Can not generate a molecular structure from the input SMILES.")
