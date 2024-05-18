from tkinter import colorchooser

def set_background_color(ax, plot_canvas):
    color = colorchooser.askcolor()[1]
    if color:
        # Set the background color of the entire plot
        fig = ax.figure
        fig.set_facecolor(color)

        # Set the background color of the canvas containing the plot
        plot_canvas.get_tk_widget().configure(bg=color)
        
        # Set the background color of each axis
        for ax in fig.get_axes():
            ax.set_facecolor(color)

        plot_canvas.draw_idle()