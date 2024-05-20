import tkinter as tk
from project_ppchem_tools_kit.set_axes_color import set_axes_color
from project_ppchem_tools_kit.set_grid_color import set_grid_color
from project_ppchem_tools_kit.set_label_color import set_label_color
from project_ppchem_tools_kit.set_background_color import set_background_color

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
