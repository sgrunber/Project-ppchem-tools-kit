import tkinter as tk
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
from tkinter import messagebox

from project_ppchem_tools_kit.set_custom_labels_and_title import set_custom_labels_and_title
from project_ppchem_tools_kit.add_data_point import add_data_point
from project_ppchem_tools_kit.set_scale import set_scale
from project_ppchem_tools_kit.toggle_grid import toggle_grid
from project_ppchem_tools_kit.display_max_point_and_coords import display_max_point_and_coords
from project_ppchem_tools_kit.open_graph_settings_window import open_graph_settings_window

graph_window = None

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
