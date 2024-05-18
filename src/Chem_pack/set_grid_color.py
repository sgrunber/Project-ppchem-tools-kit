from tkinter import colorchooser

def set_grid_color(ax, plot_canvas):
    color = colorchooser.askcolor()[1]
    if color:
        ax.grid(color=color)
        plot_canvas.draw_idle()