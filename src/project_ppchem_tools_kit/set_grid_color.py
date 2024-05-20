from tkinter import colorchooser

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