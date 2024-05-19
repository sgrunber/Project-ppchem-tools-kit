import tkinter as tk
from tkinter import Button
from project_ppchem_tools_kit.set_axes_color import set_axes_color
from project_ppchem_tools_kit.set_grid_color import set_grid_color
from project_ppchem_tools_kit.set_label_color import set_label_color
from project_ppchem_tools_kit.set_background_color import set_background_color

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
