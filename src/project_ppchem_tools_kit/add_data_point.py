from tkinter import colorchooser

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