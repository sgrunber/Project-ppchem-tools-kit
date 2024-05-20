from tkinter import colorchooser

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

