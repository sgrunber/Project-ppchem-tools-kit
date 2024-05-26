#In order to use the Tools Kit, the code must be run from its notebook file. Nonetheless, a raw code .py file containing the entirety of the code has been provided. This .py file includes all the necessary functions, classes, and logic implemented in the notebook. By reviewing the .py file, you can see the complete, unformatted version of the code in one place. This can be useful for those who prefer working directly with Python scripts or need to integrate the code into other projects. However, for the best results and to ensure all dependencies and inline outputs are correctly handled, please run the test from the provided notebook file.

import tkinter as tk
from tkinter import simpledialog, colorchooser
from sklearn.linear_model import LinearRegression
from sklearn.metrics import r2_score
from pathlib import Path
import pubchempy as pcp


from tkinter import *
import numpy as np
from tkinter import messagebox
from tkinter import filedialog
from tkinter import colorchooser

from rdkit import Chem
from rdkit.Chem import Draw, Descriptors
from rdkit.Chem.Draw.rdMolDraw2D import *
from rdkit.Chem import rdDepictor
rdDepictor.SetPreferCoordGen(True)
from rdkit import Chem

import matplotlib.pyplot as plt
import pandas as pd
from tkinter import Canvas, Entry, Button, PhotoImage
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
import matplotlib.pyplot as plt
from matplotlib_inline.backend_inline import set_matplotlib_formats
from datetime import datetime
set_matplotlib_formats('png', 'pdf')
plt.rcParams['font.family'] = 'Times New Roman'

import pyperclip

set_matplotlib_formats('png', 'pdf')
plt.rcParams['font.family'] = 'Times New Roman'

from pathlib import Path
import os
OUTPUT_PATH = Path(os.getcwd())
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
       
        df = pd.read_excel(entry_input.get().strip())
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
    """Performs a linear regression from an Excel file and plots the graph and R^2 value.

    Args:
        file_path (str): Path to the Excel file containing the data.
        x_label (str): x-axis label.
        y_label (str): y-axis label.
        title (str): Title of the graph.
        grid (bool): Whether to display grid lines on the plot. Defaults to True.
        save_as (str, optional): File path to save the graph. Defaults to None, in which case it is not saved.
        line_style (str, optional): Style of the line plot. Defaults to '-'.
        line_color (str, optional): Color of the plot line. Defaults to 'k' (black).

    Raises:
        FileNotFoundError: If the specified file_path does not exist.
    """
    def display_max_point(color):
        """Displays the maximum point on the plot with its coordinates.

        Args:
            color (str): Color of the marker and text.

        Raises:
            None
        """
        max_index = np.argmax(ax.lines[0].get_ydata())
        
        max_x = ax.lines[0].get_xdata()[max_index]
        max_y = ax.lines[0].get_ydata()[max_index]


        x_lim = ax.get_xlim()
        y_lim = ax.get_ylim()
        
        if x_lim[0] < max_x < x_lim[1] and y_lim[0] < max_y < y_lim[1]:
            x_rel = (max_x - x_lim[0]) / (x_lim[1] - x_lim[0])
            y_rel = (max_y - y_lim[0]) / (y_lim[1] - y_lim[0])
            
            dx = 0.03 * (x_lim[1] - x_lim[0])
            dy = 0.03 * (y_lim[1] - y_lim[0])
            
            if x_rel > 0.95:
                ha = 'right'
                dx = -dx
            else:
                ha = 'left'
                
            if y_rel > 0.95:
                va = 'top'
                dy = -dy
            else:
                va = 'bottom'
                
            ax.scatter(max_x, max_y, color=color, zorder=10)
            ax.text(max_x + dx, max_y + dy, f'({max_x:.1f}, {max_y:.1f})', fontsize=14, fontweight='bold', ha=ha, va=va, color=color, bbox=dict(facecolor='white', edgecolor=color, boxstyle='round,pad=0.3'))
            plot_canvas.draw_idle()
        else:
            messagebox.showwarning("Warning", "Max point coordinates are out of range.")

    def choose_color():
        """Opens a color chooser dialog and calls display_max_point with the chosen color.

        Args:
            None

        Raises:
            None
        """
        color = colorchooser.askcolor(color='red')[1]  
        if color:
            display_max_point(color)

    choose_color()


def add_data_point(ax, plot_canvas):
    """Adds data points to the graph with a chosen color.

    Args:
        ax (matplotlib.axes.Axes): The Axes object containing the plot.
        plot_canvas (FigureCanvasTkAgg): The canvas on which the plot is drawn.
    """
    color = colorchooser.askcolor()[1]
    if color:
        x_data = ax.lines[0].get_xdata()
        y_data = ax.lines[0].get_ydata()

        
        ax.scatter(x_data, y_data, color=color, zorder=10) 

        plot_canvas.draw_idle()

def set_line_color(ax, plot_canvas):
    """Sets the line color and style for the plot.

    Args:
        ax (matplotlib.axes.Axes): The Axes object containing the plot.
        plot_canvas (FigureCanvasTkAgg): The canvas on which the plot is drawn.
    """
    color = colorchooser.askcolor()[1]
    line_style = simpledialog.askstring("Line Style", "Enter line style ('-', '--', ':', '-.'):")
    if color and line_style:
        ax.lines[0].set_color(color) 
        ax.lines[0].set_linestyle(line_style)  
        plot_canvas.draw_idle()

def set_axes_color(ax, plot_canvas):
    """Sets the color of the axes and their tick marks.

    Args:
        ax (matplotlib.axes.Axes): The Axes object containing the plot.
        plot_canvas (FigureCanvasTkAgg): The canvas on which the plot is drawn.
    """
    color = colorchooser.askcolor()[1]
    if color:
        for spine in ax.spines.values():
            spine.set_color(color)
        
        ax.tick_params(axis='x', colors=color)
        ax.tick_params(axis='y', colors=color)

        plot_canvas.draw_idle()


def set_grid_color(ax, plot_canvas):
    """Sets the grid line color for the plot.

    Args:
        ax (matplotlib.axes.Axes): The Axes object containing the plot.
        plot_canvas (FigureCanvasTkAgg): The canvas on which the plot is drawn.
    """
    color = colorchooser.askcolor()[1]
    if color:
        ax.grid(color=color)
        plot_canvas.draw_idle()

def set_label_color(ax, plot_canvas):
    """Sets the color of the axis labels and title.

    Args:
        ax (matplotlib.axes.Axes): The Axes object containing the plot.
        plot_canvas (FigureCanvasTkAgg): The canvas on which the plot is drawn.
    """
    color = colorchooser.askcolor()[1]
    if color:
        ax.xaxis.label.set_color(color)
        ax.yaxis.label.set_color(color)
        ax.title.set_color(color)
        plot_canvas.draw_idle()

def set_background_color(ax, plot_canvas):
    """Sets the background color of the plot and canvas.

    Args:
        ax (matplotlib.axes.Axes): The Axes object containing the plot.
        plot_canvas (FigureCanvasTkAgg): The canvas on which the plot is drawn.
    """
    color = colorchooser.askcolor()[1]
    if color:

        fig = ax.figure
        fig.set_facecolor(color)

        # Set the background color of the canvas containing the plot
        plot_canvas.get_tk_widget().configure(bg=color)
        
        # Set the background color of each axis
        for ax in fig.get_axes():
            ax.set_facecolor(color)

        plot_canvas.draw_idle()
    
def open_graph_settings_window(fig, plot_canvas):
    """Opens a new window with options to customize the graph settings.

    Args:
        fig (matplotlib.figure.Figure): The Figure object containing the plot.
        plot_canvas (FigureCanvasTkAgg): The canvas on which the plot is drawn.
    """
    graph_settings_window = tk.Toplevel()
    graph_settings_window.title("Graph Settings")

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

    close_button = tk.Button(graph_settings_window, text="Close", command=graph_settings_window.destroy)
    close_button.pack(side=tk.BOTTOM, padx=5, pady=5)


def toggle_grid(ax, plot_canvas):
    """Toggles the visibility of grid lines on the plot.

    Args:
        ax (matplotlib.axes.Axes): The Axes object containing the plot.
        plot_canvas (FigureCanvasTkAgg): The canvas on which the plot is drawn.
    """
    ax.grid(True)
    plot_canvas.draw_idle()

def set_scale(ax, plot_canvas):
    """Sets the scale of the x and y axes.

    Args:
        ax (matplotlib.axes.Axes): The Axes object containing the plot.
        plot_canvas (FigureCanvasTkAgg): The canvas on which the plot is drawn.
    """
    x_spacing = simpledialog.askfloat("X Spacing", "Enter spacing between x ticks:")
    y_spacing = simpledialog.askfloat("Y Spacing", "Enter spacing between y ticks:")
    if x_spacing and y_spacing:
        ax.xaxis.set_major_locator(plt.MultipleLocator(x_spacing))
        ax.yaxis.set_major_locator(plt.MultipleLocator(y_spacing))
        plt.tight_layout()
        plot_canvas.draw_idle()

def set_custom_labels_and_title(ax, plot_canvas):
    """Sets custom labels and title for the plot.

    Args:
        ax (matplotlib.axes.Axes): The Axes object containing the plot.
        plot_canvas (FigureCanvasTkAgg): The canvas on which the plot is drawn.
    """
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
    """Displays the graph in a new window with options to customize it.

    Args:
        fig (matplotlib.figure.Figure): The Figure object containing the plot.
    """
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
        """Handle the closing event for the graph window.

        This function is triggered when the user attempts to close the graph window.
        It prompts the user with a confirmation dialog to ensure they want to quit.
        If the user confirms (clicks 'OK'), the graph window is hidden and the 
        main event loop is terminated.

        Args:
            None

        Returns:
            None
        """
        if messagebox.askokcancel("Quit", "Are you sure you want to quit?"):
            graph_window.withdraw()  # Hide the window
            graph_window.quit()  # Quit the main loop

            
    graph_window.protocol("WM_DELETE_WINDOW", on_closing_graph)
    graph_window.mainloop()


def calculate_error_propagation(derivatives, uncertainties):
    """
    Calculate the propagation of error using given derivatives and uncertainties.

    This function computes the standard deviation and the average value based on
    the derivatives and uncertainties provided. It follows the formula for error
    propagation in a system where each variable contributes to the overall uncertainty.

    Args:
        derivatives (list of float): The derivatives with respect to each variable.
        uncertainties (list of float): The uncertainties associated with each variable.

    Returns:
        tuple: A tuple containing the standard deviation (float) and the average value (float).
    """
    average_values = sum(derivatives) / len(derivatives)
    error = sum((unc / deriv) ** 2 for deriv, unc in zip(derivatives, uncertainties))
    standard_dev = (error ** 0.5) * average_values
    return standard_dev, average_values

def error_calculation_interface():
    """
    Provide a graphical interface for error calculation and LaTeX code generation.

    This function sets up a Tkinter window where users can input variable names,
    values, and uncertainties. It calculates the propagation of uncertainty and 
    displays the result. It also generates LaTeX code for the calculated values, 
    which can be copied to the clipboard.

    The interface includes:
        - Adding and removing variable entries.
        - Calculating and displaying the propagation of uncertainty.
        - Generating and copying LaTeX code for the results.

    Inner Functions:
        copy_latex_code(latex_code): Copy generated LaTeX code to the clipboard.
        calculate_and_display(): Calculate propagation of uncertainty and display results.
        add_variable_entry(): Add a new row of input fields for a variable.
        remove_variable_entry(row): Remove a specific row of input fields.
    """
    def copy_latex_code(latex_code):
        """
        Copy the generated LaTeX code to the clipboard.

        Args:
            latex_code (str): The LaTeX code to be copied.
        """
        pyperclip.copy(latex_code)
        result_window.clipboard_clear()
        result_window.clipboard_append(latex_code)
        result_window.update()

    def calculate_and_display():
        """
        Calculate the propagation of uncertainty and display the results.

        This function retrieves data from the input fields, calculates the error 
        propagation and average value, and displays the results in a new window. 
        It also generates the corresponding LaTeX code for the calculations.
        """
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

        error, average = calculate_error_propagation(values, uncertainties)

        result_text = f"Propagation of Uncertainty (∆z): {error:.6f}\n"
        result_text += f"Average Value: {average:.6f}\n\n"
        result_text += "Calculation Details:\n"
        for var_name, var_value, var_uncertainty in data:
            result_text += f"{var_name}: Value: {var_value}, Uncertainty: {var_uncertainty}\n"

        latex_code = r"""

        \subsection*{Propagation of Uncertainty ($\Delta\bar{z}$)}

        The propagation of uncertainty is calculated as follows:

        \begin{equation}
        \Delta \bar{z} = \bar{z} \times \sqrt{\left(\left(\frac{\Delta x}{x}\right)^2 + \left(\frac{\Delta y}{y}\right)^2 + \ldots\right)}
        \end{equation}

        The calculated uncertainty is:

        \[
        \Delta \bar{z} = """ + f"{error:.6f}" + r"""
        \]

        The average value is calculated as follows:

        \begin{equation}
        \bar{z} = \frac{1}{N} \sum_{i=1}^{N} z_i
        \end{equation}

        The calculated average value is:

        \[
        \bar{z} = """ + f"{average:.6f}" + r"""
        \]

        \subsection*{Details of Variables}

        \begin{itemize}
        """ + "\n".join([f"    \\item {var_name}: Value: {var_value}, Uncertainty: {var_uncertainty}" for var_name, var_value, var_uncertainty in data]) + r"""
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
        """
        Add a new row of input fields for a variable.

        This function dynamically adds input fields for a new variable's name, 
        value, and uncertainty to the interface.
        """
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
        """
        Remove a specific row of input fields.

        This function removes the input fields for a variable specified by the row number.

        Args:
            row (int): The row number of the variable to be removed.
        """
        entry = entries.pop(row - 1)
        for widget in entry:
            widget.destroy()
        for label in labels[(row - 1) * 3:(row - 1) * 3 + 3]:
            label.destroy()
        remove_buttons.pop(row - 1).destroy()

    error_calc_window = tk.Tk()
    error_calc_window.title("Error calculation")
    error_calc_window.configure(bg="#1B262C")

    add_variable_button = tk.Button(error_calc_window, text="Add a variable", command=add_variable_entry, height=2,
                                    bg="#0000FF", fg="#1B262C", font=("Times New Roman", 15))
    add_variable_button.grid(row=0, column=2, columnspan=3, padx=5, pady=5, sticky="")

    calculate_error_button = tk.Button(error_calc_window, text="Calcule the error propagation",
                                        command=calculate_and_display, height=2, bg="#008000", fg="#1B262C",
                                        font=("Times New Roman", 15))
    calculate_error_button.grid(row=20, column=2, columnspan=3, padx=5, pady=5, sticky="")

    # List to stock input widgets
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
        return
    elif selected_radio.get() == "5": 
        display_molecule(entry_input)
        return 

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
    """Selects all text in a widget when an event is triggered.

    Args:
        event: The event that triggered the function.

    Returns:
        str: "break" to prevent further propagation of the event.
    """
    if window and window.winfo_exists():
        event.widget.tag_add("sel", "1.0", "end")
    return "break"


def copy_text(event):
    """Copies selected text to the clipboard. Bound to an event.

    Args:
        event: The event that triggered the function.

    Returns:
        str: "break" to prevent further propagation of the event.
    """
    try:
        if window and window.winfo_exists():
            if event and event.widget and event.widget.winfo_exists():
                event.widget.event_generate("<<Copy>>")
    except Exception as e:
        print("Error while copying text:", e)
    return "break"


def bind_enter(event):
    """Binds the Enter key to the process_input function.

    Args:
        event: The event that triggered the function.
    """
    try:
        if event.keysym == "Return":
            process_input()
    except Exception as e:
        print("Error while binding Enter key:", e)

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
    """Displays a welcome message with project information in a new window."""
    global entry_input
    welcome_window = tk.Toplevel(window)
    welcome_window.title("Welcome Message")

    welcome_text = (
        "✨ Welcome to Our Project! ✨\n\n"
        "We are excited to share our project with you.\n\n"
        "Discover all the details on our GitHub repository:\n"
        "🔗 [Project ppchem Tools Kit] (https://github.com/sgrunber/Project-ppchem-tools-kit)\n\n"
        "Happy exploring! 🚀\n\n"
    )

    welcome_textbox = tk.Text(welcome_window, wrap="word", font=("Times New Roman", 20), fg="#BBE1FA", bg = "#1B262C", height=10, width=50)
    welcome_textbox.insert("1.0", welcome_text)
    welcome_textbox.config(state="disabled")
    welcome_textbox.tag_configure("center", justify="center")
    welcome_textbox.tag_add("center", "1.0", "end")
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
    """Updates the clock label with the current time and date every second."""
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
    """
    Handle the closing event for the graph window.

    This function is triggered when the user attempts to close the graph window. 
    It prompts the user with a confirmation dialog to ensure they want to quit. 
    If the user confirms (clicks 'OK'), the graph window is hidden and the 
    main event loop is terminated.
    """
    global window
    if messagebox.askokcancel("Quit", "Are you sure you want to quit?"):
        window.destroy()
window.protocol("WM_DELETE_WINDOW", on_closing)

window.resizable(False, False)

window.bind("<Return>", process_input)

window.mainloop()
