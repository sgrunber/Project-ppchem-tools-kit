from tkinter import colorchooser

def set_label_color(ax, plot_canvas):
    color = colorchooser.askcolor()[1]
    if color:
        ax.xaxis.label.set_color(color)
        ax.yaxis.label.set_color(color)
        ax.title.set_color(color)
        plot_canvas.draw_idle()