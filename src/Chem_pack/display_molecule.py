from rdkit import Chem
from rdkit.Chem import Descriptors, Draw
from svglib.svglib import svg2rlg
from reportlab.graphics import renderPDF, renderPM
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk, FigureCanvasTk
from io import BytesIO
import tkinter as tk
from PIL import Image

def display_molecule(entry_input):
    
    smiles = entry_input.get().strip()  # Obtenir la chaîne SMILES à partir de l'entrée utilisateur
    mol = Chem.MolFromSmiles(smiles)  # Générer la structure moléculaire à partir de la chaîne SMILES

    if mol is not None:
        # Fonction pour afficher la molécule en 2D
        def show_molecule_2d():
            def get_svg(mol):
                d2d = Draw.MolDraw2DSVG(350, 300)
                d2d.DrawMolecule(mol)
                svg_data = d2d.GetDrawingText()
                drawing = svg2rlg(BytesIO(svg_data.encode()))
                renderPM.drawToFile(drawing, "temp.png", fmt="PNG")
            get_svg(mol)
            img = Image.open('temp.png')

            # Créer une nouvelle figure matplotlib
            fig, ax = plt.subplots(figsize=(4, 3))  # Définir la taille de la figure
            ax.imshow(img)  # Afficher l'image sur la figure
            ax.axis('off')  # Désactiver les axes
            plt.tight_layout()
            # Créer une fenêtre Tkinter pour afficher la structure moléculaire en 2D
            
            mol_window = tk.Toplevel()
            mol_window.title("Molecular Structure")

            # Convertir la figure en format Tkinter PhotoImage
            pimg = FigureCanvasTkAgg(fig, master=mol_window)
            pimg.draw()

            # Créer un canevas Tkinter pour afficher l'image
            canvas = pimg.get_tk_widget()
            canvas.pack()

            # Ajouter la barre d'outils de navigation
            toolbar = NavigationToolbar2Tk(pimg, mol_window)  # Passer la figure à la barre d'outils
            toolbar.update()
            toolbar.pack()
            # Enregistrer l'image avec une meilleure qualité


            # Ajouter une barre d'outils pour enregistrer l'image


            mol_window.mainloop()

        show_molecule_2d()
    else:
        print("Erreur : Impossible de générer une structure moléculaire à partir du SMILES fourni.")
