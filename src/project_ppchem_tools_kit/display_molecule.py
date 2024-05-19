from rdkit import Chem
from rdkit.Chem import Draw
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
import tkinter as tk



def display_molecule(entry_input):
    
    smiles = entry_input.get().strip()  # Obtenir la chaîne SMILES à partir de l'entrée utilisateur
    mol = Chem.MolFromSmiles(smiles)  # Générer la structure moléculaire à partir de la chaîne SMILES

    if mol is not None:
        # Fonction pour afficher la molécule en 2D
        def show_molecule_2d():
            img = Draw.MolToImage(mol)

            # Créer une nouvelle fenêtre Tkinter pour afficher la structure moléculaire en 2D
            mol_window = tk.Toplevel()
            mol_window.title("Molecular Structure")

            # Convertir l'image en format Tkinter PhotoImage
            pimg = FigureCanvasTkAgg(plt.Figure(figsize=(4, 3)), master=mol_window)
            ax = pimg.figure.add_subplot(111)
            ax.imshow(img, interpolation='bilinear')
            ax.axis('off')

            # Créer un canevas Tkinter pour afficher l'image
            canvas = pimg.get_tk_widget()
            canvas.pack()

            # Ajouter la barre d'outils de navigation
            toolbar = NavigationToolbar2Tk(pimg, mol_window)  # Passer la figure à la barre d'outils
            toolbar.update()
            toolbar.pack()

            mol_window.mainloop()

        show_molecule_2d()
    else:
        print("Erreur : Impossible de générer une structure moléculaire à partir du SMILES fourni.")
