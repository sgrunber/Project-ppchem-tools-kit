import matplotlib.pyplot as plt
import pandas as pd
from tkinter import messagebox
from project_ppchem_tools_kit.graph_function import display_graph
from sklearn.linear_model import LinearRegression
from sklearn.metrics import r2_score
import numpy as np


def linear_regression(entry_input, x_label, y_label, title, grid=True, save_as=None, line_style='-', line_color='k'):
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
    try:
        df = pd.read_excel(entry_input.get().strip())  
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