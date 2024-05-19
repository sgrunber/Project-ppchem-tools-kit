from tkinter import simpledialog
import matplotlib.pyplot as plt

def set_custom_labels_and_title(ax, plot_canvas):
    x_label = simpledialog.askstring("Custom Labels and Title", "Enter X-axis Label:")
    y_label = simpledialog.askstring("Custom Labels and Title", "Enter Y-axis Label:")
    title = simpledialog.askstring("Custom Labels and Title", "Enter Graph Title:")
    if x_label and y_label and title:
        ax.set_xlabel(x_label, fontsize=20)
        ax.set_ylabel(y_label, fontsize=20)
        ax.set_title(title, fontsize=25, fontweight='bold')
        plt.tight_layout() 
        plot_canvas.draw_idle()