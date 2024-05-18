from tkinter import colorchooser

def set_axes_color(ax, plot_canvas):
    color = colorchooser.askcolor()[1]
    if color:
        for spine in ax.spines.values():
            spine.set_color(color)
        
        # Changer la couleur des valeurs sur les axes
        ax.tick_params(axis='x', colors=color)
        ax.tick_params(axis='y', colors=color)

        plot_canvas.draw_idle()

