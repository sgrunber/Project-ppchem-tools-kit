import numpy as np
import matplotlib.pyplot as plt
import tkinter as tk
from tkinter import messagebox
from tkinter import simpledialog, colorchooser
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk

#NOTE : In order for the code to run, graph_function.py englobes the following functions :
#   display_max_point_and_coords, 
#   add_data_point, 
#   set_line_color, 
#   set_axes_color, 
#   set_grud_color, 
#   set_label_color, 
#   set_background_color, 
#   open_graph_settings_window, 
#   toggle_grid, 
#   set_scale, 
#   set_custom_labels_and_title, 
#   display_graph, 
#   on_closing_graph
# as well as the nested functions. However, in order to test each of the functions listed hereabove, separate _.py files were created for each main function.

graph_window = None 

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
            
            dx = 0.01 * (x_lim[1] - x_lim[0])
            dy = 0.05 * (y_lim[1] - y_lim[0])
            
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

        plot_canvas.get_tk_widget().configure(bg=color)
        
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

    set_line_color_button = tk.Button(graph_settings_window, text="Set Line Color", command=lambda: set_line_color(fig.axes[0], plot_canvas))
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
        ax.lines[0].set_linewidth(1) 


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

