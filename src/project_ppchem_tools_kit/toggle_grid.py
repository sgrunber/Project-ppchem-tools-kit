def toggle_grid(ax, plot_canvas):
    """Toggles the visibility of grid lines on the plot.

    Args:
        ax (matplotlib.axes.Axes): The Axes object containing the plot.
        plot_canvas (FigureCanvasTkAgg): The canvas on which the plot is drawn.
    """
    ax.grid(True)
    plot_canvas.draw_idle()
