
import matplotlib.pyplot as plt
import pandas as pd
from tkinter import messagebox
from project_ppchem_tools_kit.graph_function import display_graph

def make_graph(entry_input, x_label, y_label, title, grid=True, save_as=None, line_style='-', line_color='k'):
    try:
        # Récupérez les données du fichier Excel ou d'une autre source
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
        
