
import matplotlib.pyplot as plt
import pandas as pd
from tkinter import messagebox
from project_ppchem_tools_kit.graph_function import display_graph

def make_graph(entry_input, x_label, y_label, title, grid=True, save_as=None, line_style='-', line_color='k'):
    """Creates a scatter plot from data in a given Excel file.

    Args:
        filepath (str): Path to the Excel file containing the data.
        x_label (str): x-axis label.
        y_label (str): y-axis label.
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
        
