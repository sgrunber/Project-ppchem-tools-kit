from tkinter import colorchooser, simpledialog

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