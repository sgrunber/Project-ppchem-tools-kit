from tkinter import colorchooser, simpledialog

def set_line_color(ax, plot_canvas):
    color = colorchooser.askcolor()[1]
    line_style = simpledialog.askstring("Line Style", "Enter line style ('-', '--', ':', '-.'):")
    if color and line_style:
        ax.lines[0].set_color(color)  # Change seulement la couleur de la première ligne (la courbe)
        ax.lines[0].set_linestyle(line_style)  
        plot_canvas.draw_idle()