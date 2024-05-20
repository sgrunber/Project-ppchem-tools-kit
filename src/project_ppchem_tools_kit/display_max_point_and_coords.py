import numpy as np
from tkinter import messagebox, colorchooser


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