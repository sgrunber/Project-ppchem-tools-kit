

from tkinter import ttk
import tkinter as tk
from tkinter import simpledialog, colorchooser
from sklearn.linear_model import LinearRegression
from sklearn.metrics import r2_score
from pathlib import Path
from IPython.display import display, Image
import os
import pubchempy as pcp


from tkinterweb import HtmlFrame

from tkinter import *
from matplotlib import colors as mcolors
import numpy as np
from PIL import Image, ImageDraw, ImageTk
from tkinter import messagebox
from tkinter import filedialog
from tkinter import colorchooser
import subprocess

from rdkit import Chem
from rdkit.Chem import Draw
import tempfile
import requests
from rdkit.Chem import ChemicalFeatures, MolFromSmiles, Draw
from rdkit.Chem import Draw, Lipinski, Crippen, Descriptors
from rdkit.Chem import ChemicalFeatures
from rdkit.Chem.Draw import rdMolDraw2D
from rdkit.Chem.Draw.rdMolDraw2D import *
from rdkit.Chem import rdDepictor
rdDepictor.SetPreferCoordGen(True)
from rdkit.Chem.Draw import IPythonConsole
from IPython.display import SVG
from rdkit import Chem
from rdkit.Chem import Draw, AllChem
from rdkit.Chem.Draw import rdMolDraw2D
from collections import defaultdict

import matplotlib.pyplot as plt
import pandas as pd
from tkinter import Tk, Canvas, Entry, Text, Button, PhotoImage
from tkinter import IntVar, Radiobutton
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk, FigureCanvasTk
import matplotlib.pyplot as plt
from matplotlib_inline.backend_inline import set_matplotlib_formats
from matplotlib.figure import Figure
from datetime import datetime
set_matplotlib_formats('png', 'pdf')
plt.rcParams['font.family'] = 'Times New Roman'
import sys
sys.path.insert(0, "../src")
#from "src/name_to_smiles" import name_to_smiles``

import pyperclip

'''from Chem_pack.smiles_to_molar_mass import smiles_to_molar_mass
from Chem_pack.name_to_smiles import name_to_smiles
from Chem_pack.display_molecule import display_molecule'''

set_matplotlib_formats('png', 'pdf')
plt.rcParams['font.family'] = 'Times New Roman'

OUTPUT_PATH = Path(__file__).parent
ASSETS_PATH = OUTPUT_PATH / Path("./assets/frame0")

def relative_to_assets(path: str) -> Path:
    """Constructs a path to a file located in the assets directory by combining the provided relative path with the ASSETS_PATH
    
    Args : 
        path (str) ; the relative path to be appended
    
    Returns : 
        the full path to the asset
    """
    return ASSETS_PATH / Path(path)

##############################################################################################################################################################

entry_input = None
selected_radio = None
window = None
graph_window = None
mol_window = None
x_label = None
y_label = None
title = None
result_text = None


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
       
        df = pd.read_excel(entry_input.get().strip())  # Utilisez entry_input pour obtenir le chemin du fichier
        data_list = df.values.tolist()
        x_values = [row[0] for row in data_list]
        y_values = [row[1] for row in data_list]

        fig, ax = plt.subplots()
        plt.plot(x_values, y_values, linestyle=line_style, color=line_color)
        plt.xlabel(x_label, fontsize=20)
        plt.ylabel(y_label, fontsize=20)
        plt.title(title, fontsize=25, fontweight='bold')
        plt.tight_layout()

        if save_as:
            plt.savefig(save_as)

        display_graph(fig)

    except FileNotFoundError:
        messagebox.showerror("Error", "The data file was not found.")



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
        # Récupérez les données du fichier Excel ou d'une autre source
        df = pd.read_excel(entry_input.get().strip())  # Utilisez entry_input pour obtenir le chemin du fichier
        data_list = df.values.tolist()
        x_values = np.array([row[0] for row in data_list]).reshape(-1, 1)
        y_values = [row[1] for row in data_list]

        model = LinearRegression()
        model.fit(x_values, y_values)

        y_pred = model.predict(x_values)
        r2 = r2_score(y_values, y_pred)

        fig, ax = plt.subplots()
        plt.plot(x_values, y_pred, color='red', label='Linear Regression', linestyle=line_style, linewidth=1)
        plt.scatter(x_values, y_values, color='blue', label='Data Points') # peut-être enlevé les label
        plt.xlabel(x_label, fontsize=20)
        plt.ylabel(y_label, fontsize=20)
        plt.title(title, fontsize=25, fontweight='bold')
        plt.tight_layout()
        plt.legend()
        plt.text(0.6, 0.8, f'$R^2 = {r2:.2f}$', ha='center', va='center', transform=ax.transAxes, fontsize=13, fontname='Times New Roman')

        if save_as:
            plt.savefig(save_as)

        display_graph(fig)

    except FileNotFoundError:
        messagebox.showerror("Error", "The data file was not found.")

def display_max_point_and_coords(ax, plot_canvas):
    def display_max_point(color):
        # Trouver l'indice du maximum de la courbe
        max_index = np.argmax(ax.lines[0].get_ydata())
        
        # Récupérer les coordonnées du maximum
        max_x = ax.lines[0].get_xdata()[max_index]
        max_y = ax.lines[0].get_ydata()[max_index]

        # Déterminer la position relative du maximum par rapport aux axes
        x_lim = ax.get_xlim()
        y_lim = ax.get_ylim()
        
        if x_lim[0] < max_x < x_lim[1] and y_lim[0] < max_y < y_lim[1]:
            # Déterminer la position relative du point sur les axes
            x_rel = (max_x - x_lim[0]) / (x_lim[1] - x_lim[0])
            y_rel = (max_y - y_lim[0]) / (y_lim[1] - y_lim[0])
            
            # Calculer les décalages en fonction de la position relative
            dx = 0.03 * (x_lim[1] - x_lim[0])
            dy = 0.03 * (y_lim[1] - y_lim[0])
            
            # Si le point est proche du bord droit, décaler à gauche
            if x_rel > 0.95:
                ha = 'right'
                dx = -dx
            else:
                ha = 'left'
                
            # Si le point est proche du bord supérieur, décaler vers le bas
            if y_rel > 0.95:
                va = 'top'
                dy = -dy
            else:
                va = 'bottom'
                
            # Afficher les coordonnées dans le graphique avec une seule décimale
            ax.scatter(max_x, max_y, color=color, zorder=10)
            ax.text(max_x + dx, max_y + dy, f'({max_x:.1f}, {max_y:.1f})', fontsize=14, fontweight='bold', ha=ha, va=va, color=color, bbox=dict(facecolor='white', edgecolor=color, boxstyle='round,pad=0.3'))
            plot_canvas.draw_idle()
        else:
            messagebox.showwarning("Warning", "Max point coordinates are out of range.")

    def choose_color():
        color = colorchooser.askcolor(color='red')[1]  # Choisissez la couleur par défaut ici
        if color:
            display_max_point(color)

    choose_color()

def add_data_point(ax, plot_canvas):
    # Demander à l'utilisateur de choisir une couleur
    color = colorchooser.askcolor()[1]
    if color:
        # Récupérer les données x et y du graphique
        x_data = ax.lines[0].get_xdata()
        y_data = ax.lines[0].get_ydata()

        # Ajouter les points avec la couleur spécifiée à l'aide de ax.scatter
        
        ax.scatter(x_data, y_data, color=color, zorder=10) 

        # Mettre à jour le canvas pour afficher les modifications
        plot_canvas.draw_idle()

def set_line_color(ax, plot_canvas):
    color = colorchooser.askcolor()[1]
    line_style = simpledialog.askstring("Line Style", "Enter line style ('-', '--', ':', '-.'):")
    if color and line_style:
        ax.lines[0].set_color(color)  # Change seulement la couleur de la première ligne (la courbe)
        ax.lines[0].set_linestyle(line_style)  
        plot_canvas.draw_idle()


def set_axes_color(ax, plot_canvas):
    color = colorchooser.askcolor()[1]
    if color:
        for spine in ax.spines.values():
            spine.set_color(color)
        
        # Changer la couleur des valeurs sur les axes
        ax.tick_params(axis='x', colors=color)
        ax.tick_params(axis='y', colors=color)

        plot_canvas.draw_idle()


def set_grid_color(ax, plot_canvas):
    color = colorchooser.askcolor()[1]
    if color:
        ax.grid(color=color)
        plot_canvas.draw_idle()

def set_label_color(ax, plot_canvas):
    color = colorchooser.askcolor()[1]
    if color:
        ax.xaxis.label.set_color(color)
        ax.yaxis.label.set_color(color)
        ax.title.set_color(color)
        plot_canvas.draw_idle()

def set_background_color(ax, plot_canvas):
    color = colorchooser.askcolor()[1]
    if color:
        # Set the background color of the entire plot
        fig = ax.figure
        fig.set_facecolor(color)

        # Set the background color of the canvas containing the plot
        plot_canvas.get_tk_widget().configure(bg=color)
        
        # Set the background color of each axis
        for ax in fig.get_axes():
            ax.set_facecolor(color)

        plot_canvas.draw_idle()
    
def open_graph_settings_window(fig, plot_canvas):
    graph_settings_window = tk.Toplevel()
    graph_settings_window.title("Graph Settings")

    # Boutons pour les paramètres du graphique
    set_line_color_button = tk.Button(graph_settings_window, text="Set Line Color", command=lambda: set_axes_color(fig.axes[0], plot_canvas))
    set_line_color_button.pack(side=tk.TOP, padx=5, pady=5)

    set_axes_color_button = tk.Button(graph_settings_window, text="Set Axis Color", command=lambda: set_axes_color(fig.axes[0], plot_canvas))
    set_axes_color_button.pack(side=tk.TOP, padx=5, pady=5)

    set_grid_color_button = tk.Button(graph_settings_window, text="Set grid", command=lambda: set_grid_color(fig.axes[0], plot_canvas))
    set_grid_color_button.pack(side=tk.TOP, padx=5, pady=5)

    set_label_color_button = tk.Button(graph_settings_window, text="Set Label Color", command=lambda: set_label_color(fig.axes[0], plot_canvas))
    set_label_color_button.pack(side=tk.TOP, padx=5, pady=5)

    set_background_color_button = tk.Button(graph_settings_window, text="Set Background Color", command=lambda: set_background_color(fig.axes[0], plot_canvas))
    set_background_color_button.pack(side=tk.TOP, padx=5, pady=5)

    # Bouton pour fermer la fenêtre des paramètres du graphique
    close_button = tk.Button(graph_settings_window, text="Close", command=graph_settings_window.destroy)
    close_button.pack(side=tk.BOTTOM, padx=5, pady=5)


def toggle_grid(ax, plot_canvas):
    # Set the grid lines visibility to True
    ax.grid(True)
    plot_canvas.draw_idle()

def set_scale(ax, plot_canvas):
    x_spacing = simpledialog.askfloat("X Spacing", "Enter spacing between x ticks:")
    y_spacing = simpledialog.askfloat("Y Spacing", "Enter spacing between y ticks:")
    if x_spacing and y_spacing:
        ax.xaxis.set_major_locator(plt.MultipleLocator(x_spacing))
        ax.yaxis.set_major_locator(plt.MultipleLocator(y_spacing))
        plt.tight_layout()
        plot_canvas.draw_idle()

def set_custom_labels_and_title(ax, plot_canvas):
    x_label = simpledialog.askstring("Custom Labels and Title", "Enter X-axis Label:")
    y_label = simpledialog.askstring("Custom Labels and Title", "Enter Y-axis Label:")
    title = simpledialog.askstring("Custom Labels and Title", "Enter Graph Title:")
    if x_label and y_label and title:
        ax.set_xlabel(x_label, fontsize=20)
        ax.set_ylabel(y_label, fontsize=20)
        ax.set_title(title, fontsize=25, fontweight='bold')
        plt.tight_layout() 
        plot_canvas.draw_idle()
        


def display_graph(fig):
    global graph_window
    if graph_window:
        graph_window.destroy()

    graph_window = tk.Toplevel()
    graph_window.title("Graph")
    graph_window.geometry("1000x600")


    plot_canvas = FigureCanvasTkAgg(fig, master=graph_window)
    plot_canvas.draw()
    plot_canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=1)

    toolbar = NavigationToolbar2Tk(plot_canvas, graph_window)
    toolbar.update()

    
    for ax in fig.get_axes():
        ax.set_autoscale_on(False)
        ax.set_xlim(auto=True)
        ax.set_ylim(auto=True)
        ax.lines[0].set_linewidth(1)  # Reset line width


    custom_button = tk.Button(graph_window, text="Custom Labels and Title", command=lambda: set_custom_labels_and_title(fig.axes[0], plot_canvas))
    custom_button.pack(side=tk.LEFT, padx=5)

    add_data_point_button = tk.Button(graph_window, text="Add Data Point", command=lambda: add_data_point(fig.axes[0], plot_canvas))
    add_data_point_button.pack(side=tk.LEFT, padx=5)

    set_scale_button = tk.Button(graph_window, text="Set Scale", command=lambda: set_scale(fig.axes[0], plot_canvas))
    set_scale_button.pack(side=tk.LEFT, padx=5)

    grid_button = tk.Button(graph_window, text="Toggle Grid", command=lambda: toggle_grid(fig.axes[0], plot_canvas))
    grid_button.pack(side=tk.LEFT, padx=5)
    
    display_max_button = tk.Button(graph_window, text="Display Max Point", command=lambda: display_max_point_and_coords(fig.axes[0], plot_canvas))
    display_max_button.pack(side=tk.LEFT, padx=5)

    graph_settings_button = tk.Button(graph_window, text="Graph Settings", command=lambda: open_graph_settings_window(fig, plot_canvas))
    graph_settings_button.pack(side=tk.LEFT, padx=5, pady=5)
    
    
    def on_closing_graph():
        if messagebox.askokcancel("Quit", "Are you sure you want to quit?"):
            graph_window.withdraw()  # Hide the window
            graph_window.quit()  # Quit the main loop

            
    graph_window.protocol("WM_DELETE_WINDOW", on_closing_graph)
    graph_window.mainloop()


def error_calculation_interface():
    # Fonction pour calculer la propagation d'erreur
    def calculate_error_propagation(derivatives, uncertainties):
        average_values = sum(derivatives) / len(derivatives)
        error = sum((unc / deriv) ** 2 for deriv, unc in zip(derivatives, uncertainties))
        standard_dev = (error** 0.5) * average_values
        return standard_dev, average_values 

    def copy_latex_code(latex_code):
        pyperclip.copy(latex_code)
        result_window.clipboard_clear()
        result_window.clipboard_append(latex_code)
        result_window.update()

    def calculate_and_display():
        data = []
        values = []
        uncertainties = []
        for i in range(len(entries)):
            var_name = entry_vars[i].get()
            var_value = float(entry_values[i].get())
            var_uncertainty = float(entry_uncertainties[i].get())
            data.append((var_name, var_value, var_uncertainty))
            values.append(var_value)
            uncertainties.append(var_uncertainty)

        # Calculate the propagation of uncertainty
        error, average = calculate_error_propagation(values, uncertainties)

        # Display the result with calculation details
        result_text = f"Propagation of Uncertainty (∆z): {error:.6f}\n"
        result_text += f"Average Value: {average:.6f}\n\n"
        result_text += "Calculation Details:\n"
        for var_name, var_value, var_uncertainty in data:
            result_text += f"{var_name}: Value: {var_value}, Uncertainty: {var_uncertainty}\n"

      
        latex_code = r"""


        Propagation of Uncertainty ($\Delta z$)
        The propagation of uncertainty is calculated as follows:
        \begin{equation}
        \Delta \bar{z} = \bar{z} \times \sqrt{\left(\left(\frac{\Delta x}{x}\right)^2 + \left(\frac{\Delta y}{y}\right)^2 + \ldots\right)}
        \end{equation}
        \[
        \Delta z = """ + f"{error:.6f}" + r"""
        \]
        
        The average value is calculated as follows:
        
        \begin{equation}
        \bar{z} = \frac{1}{N} \sum_{i=1}^{N} z_i
        \end{equation}
        \[
        \bar{z} = """ + f"{average:.6f}" + r"""
        \]

        The details of variables:
        \begin{itemize}
        """ + "\n".join([f"    \item {var_name}: Value: {var_value}, Uncertainty: {var_uncertainty}" for var_name, var_value, var_uncertainty in data]) + r"""
        \end{itemize}
        """

        # Create a new window to display the result
        global result_window
        result_window = tk.Toplevel(error_calc_window)
        result_window.title("Result")

        # Display the result in a Text widget
        result_textbox = tk.Text(result_window, wrap="word", font=("Times New Roman", 20), height=15, width=70)
        result_textbox.insert("1.0", result_text)
        result_textbox.grid(row=0, column=0, padx=5, pady=5)

        # Button to copy LaTeX code to clipboard
        copy_button = tk.Button(result_window, text="Copy LaTeX Code", command=lambda: copy_latex_code(latex_code))
        copy_button.grid(row=1, column=0, padx=5, pady=5)

        
    def add_variable_entry():
        # Ajouter des champs d'entrée pour une nouvelle variable
        row = len(entries) + 1

        var_name_label = tk.Label(error_calc_window, text="Variable Name", bg="#1B262C", fg="white")
        var_name_label.grid(row=row, column=0, padx=5, pady=5)
        var_name_entry = tk.Entry(error_calc_window)
        var_name_entry.grid(row=row, column=1, padx=5, pady=5)
        entry_vars.append(var_name_entry)
        labels.append(var_name_label)

        var_value_label = tk.Label(error_calc_window, text="Value", bg="#1B262C", fg="white")
        var_value_label.grid(row=row, column=2, padx=5, pady=5)
        var_value_entry = tk.Entry(error_calc_window)
        var_value_entry.grid(row=row, column=3, padx=5, pady=5)
        entry_values.append(var_value_entry)
        labels.append(var_value_label)

        var_uncertainty_label = tk.Label(error_calc_window, text="Uncertainty", bg="#1B262C", fg="white")
        var_uncertainty_label.grid(row=row, column=4, padx=5, pady=5)
        var_uncertainty_entry = tk.Entry(error_calc_window)
        var_uncertainty_entry.grid(row=row, column=5, padx=5, pady=5)
        entry_uncertainties.append(var_uncertainty_entry)
        labels.append(var_uncertainty_label)

        entries.append((var_name_entry, var_value_entry, var_uncertainty_entry))

        remove_variable_button = tk.Button(error_calc_window, text="Remove", command=lambda r=row: remove_variable_entry(r),
                                            bg="#FF0000", fg="#1B262C")
        remove_variable_button.grid(row=row, column=6, padx=5, pady=5)
        remove_buttons.append(remove_variable_button)

    def remove_variable_entry(row):
        entry = entries.pop(row - 1)
        for widget in entry:
            widget.destroy()
        for label in labels[(row - 1) * 3:(row - 1) * 3 + 3]:
            label.destroy()
        remove_buttons.pop(row - 1).destroy()

    error_calc_window = tk.Tk()
    error_calc_window.title("Calcul d'erreur")
    error_calc_window.configure(bg="#1B262C")

    add_variable_button = tk.Button(error_calc_window, text="Ajouter une variable", command=add_variable_entry, height=2,
                                    bg="#0000FF", fg="#1B262C", font=("Times New Roman", 15))
    add_variable_button.grid(row=0, column=2, columnspan=3, padx=5, pady=5, sticky="")

    # Pour le bouton "Calculer la propagation d'erreur"
    calculate_error_button = tk.Button(error_calc_window, text="Calculer la propagation d'erreur",
                                        command=calculate_and_display, height=2, bg="#008000", fg="#1B262C",
                                        font=("Times New Roman", 15))
    calculate_error_button.grid(row=20, column=2, columnspan=3, padx=5, pady=5, sticky="")

    # Listes pour stocker les widgets d'entrée
    entries = []
    entry_vars = []
    entry_values = []
    entry_uncertainties = []
    remove_buttons = []
    labels = []

    error_calc_window.mainloop()


def process_input(event=None):
    """Processes the input based on the selected mode (molecule name, SMILES, file path for plotting). Validates the input and displays the appropriate results or errors.

    Args:
        event (Event, optional): The event that triggered this function, such as a button click or another event from a GUI component. Defaults to None.

    Returns:
        None: Does not return any values, directly affects the GUI by displaying messages or updating GUI components.
    """
    global selected_radio, x_label, y_label, title, result_text
    if selected_radio.get() == "1":
        smiles_molecule = name_to_smiles(entry_input.get().strip())
        if smiles_molecule is not None:
            result_text = f"The SMILES for the molecule '{entry_input.get().strip()}' is: {smiles_molecule}"
        else:
            result_text = f"The molecule '{entry_input.get().strip()}' was not found in the PubChem database."
    elif selected_radio.get() == "2":
        molar_mass = smiles_to_molar_mass(entry_input.get().strip())
        if molar_mass is not None:
            result_text = f"The molar mass of the molecule is: {molar_mass} g/mol"
        else:
            result_text = "Invalid SMILES or molecule not found."
    elif selected_radio.get() in ["3", "4"]:
        if selected_radio.get() == "3":
            make_graph(entry_input, x_label, y_label, title)
        else:
            linear_regression(entry_input, x_label, y_label, title)
        return  # Return to prevent displaying the result window
    elif selected_radio.get() == "5":  # Ajouter cette condition pour le bouton radio "5"
        display_molecule(entry_input)  # Appeler la fonction pour afficher la structure moléculaire
        return  # Return to prevent displaying the result window

    elif selected_radio.get() == "6":
        error_calculation_interface()
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

def select_all(event):
    if window and window.winfo_exists():
        event.widget.tag_add("sel", "1.0", "end")
    return "break"


def copy_text(event):
    """Copies selected text to the clipboard. Bound to an event.

    Args:
        event: The event that triggered the function.

    Returns:
        None: Does not return any values.
    """
    try:
        if window and window.winfo_exists():
            if event and event.widget and event.widget.winfo_exists():
                event.widget.event_generate("<<Copy>>")
    except Exception as e:
        print("Error while copying text:", e)
    return "break"


# Avant de détruire la fenêtre principale, détacher les gestionnaires d'événements


def bind_enter(event):
    try:
        if event.keysym == "Return":
            process_input()
    except Exception as e:
        print("Error while binding Enter key:", e)

# Attachez le gestionnaire d'événements à la fenêtre principale
if window:
    window.bind("<KeyPress>", bind_enter)


############################################################################################################################################

window = tk.Tk()  #creates a Tkinter window instance

x_axis_label = tk.StringVar()
y_axis_label = tk.StringVar()
graph_title = tk.StringVar()

window.title("Project Tools-Kit")
input_type = tk.IntVar() #creates a Tkinter IntVar variable, which is used to track the value of the selected input type. In this case, it's initialized to 1
input_type.set(1)



def welcome_message():
    global entry_input
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

canvas = Canvas(window, bg="#1B262C", height=800, width=1000, bd=0, highlightthickness=0, relief="ridge")
canvas.place(x=0, y=0)

canvas.create_text(12.0, 194.0, anchor="nw", text="Choose Input Type :", fill="#BBE1FA", font=("Times New Roman", 27, "bold"))


canvas.create_text(92.0, 378.0, anchor="nw",text="Input :", fill="#BBE1FA", font=("Times New Roman", 30, "bold"))
entry_input = Entry(bd=1, bg="#0F4C75", fg="#BBE1FA", highlightthickness=0, font=("Times New Roman", 25))
entry_input.place(x=272.5, y=368.0, width=605.0, height=74.0)

button_image_browse = PhotoImage(file=relative_to_assets("button_browse.png")) 
button_browse = Button(image=button_image_browse, borderwidth=0, highlightthickness=0, command=browse_excel_file, relief="flat")
button_browse.place(x=409.0, y=500.0, width=154.2904052734375, height=60.0)

button_image_process = PhotoImage(file=relative_to_assets("button_process.png"))
button_process = Button(image=button_image_process, borderwidth=0, highlightthickness=0, command=process_input, relief="flat")
button_process.place(x=394.5, y=620.0, width=182.0, height=73.0)

button_image_title = PhotoImage(file=relative_to_assets("button_title.png"))
button_title = Button(image=button_image_title, borderwidth=0, highlightthickness=0, command=welcome_message, relief="flat")
button_title.place(x=273.0, y=36.0, width=455.0, height=90.0)



logo_image = PhotoImage(file=relative_to_assets('EPFL_logo.png'))
logo_label = tk.Label(window, image=logo_image)
logo_label.place(x=10.0, y=10.0, width=150.0, height=60.0)

clock_label = tk.Label(window, text="", font=("Times New Roman", 20), bg="#1B262C", fg="#BBE1FA")
clock_label.place(x=900, y=10)



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
    (291, 252, "Excel Graph", "3"),
    (700, 194, "Linear Regression", "4"), 
    (700, 252, "Show Molecule", "5"),
    (469, 252, "Error calculation", "6"),
    (469, 194, "Molecular weight", "2")
]   
selected_radio = tk.StringVar(value="none")

for data in radio_button_data:
    create_radio_button(*data)



def update_clock():
    current_time = datetime.now().strftime('%H:%M:%S')
    current_date = datetime.now().strftime('%Y-%m-%d')
    clock_label.config(text=current_time)
    date_label.config(text=current_date)
    window.after(1000, update_clock)

date_label = tk.Label(window, text="", font=("Times New Roman", 20), bg="#1B262C", fg="#FFFFFF")
date_label.place(x=890, y=10)

clock_label = tk.Label(window, text="", font=("Times New Roman", 20), bg="#1B262C", fg="#FFFFFF")
clock_label.place(x=890, y=40)

update_clock()

window.unbind_all("<Return>")

def on_closing():
    global window
    if messagebox.askokcancel("Quit", "Are you sure you want to quit?"):
        window.destroy()
window.protocol("WM_DELETE_WINDOW", on_closing)

window.resizable(False, False)

window.bind("<Return>", process_input)

window.mainloop()
