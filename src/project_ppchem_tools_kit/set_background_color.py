from tkinter import colorchooser

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

        # Set the background color of the canvas containing the plot
        plot_canvas.get_tk_widget().configure(bg=color)
        
        # Set the background color of each axis
        for ax in fig.get_axes():
            ax.set_facecolor(color)

        plot_canvas.draw_idle()