import numpy as np
import matplotlib.pyplot as plt
import tkinter as tk
from tkinter import messagebox
from tkinter import simpledialog, colorchooser
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk, FigureCanvasTk
from matplotlib.figure import Figure


graph_window = None 

def display_max_point_and_coords(ax, plot_canvas):
    def display_max_point(color):
        # Trouver l'indice du maximum de la courbe
        max_index = np.argmax(ax.lines[0].get_ydata())
        
        # Récupérer les coordonnées du maximum
        max_x = ax.lines[0].get_xdata()[max_index]
        max_y = ax.lines[0].get_ydata()[max_index]

        # Déterminer la position relative du maximum par rapport aux axes
        x_lim = ax.get_xlim()
        y_lim = ax.get_ylim()
        
        if x_lim[0] < max_x < x_lim[1] and y_lim[0] < max_y < y_lim[1]:
            # Déterminer la position relative du point sur les axes
            x_rel = (max_x - x_lim[0]) / (x_lim[1] - x_lim[0])
            y_rel = (max_y - y_lim[0]) / (y_lim[1] - y_lim[0])
            
            # Calculer les décalages en fonction de la position relative
            dx = 0.01 * (x_lim[1] - x_lim[0])
            dy = 0.05 * (y_lim[1] - y_lim[0])
            
            # Si le point est proche du bord droit, décaler à gauche
            if x_rel > 0.95:
                ha = 'right'
                dx = -dx
            else:
                ha = 'left'
                
            # Si le point est proche du bord supérieur, décaler vers le bas
            if y_rel > 0.95:
                va = 'top'
                dy = -dy
            else:
                va = 'bottom'
                
            # Afficher les coordonnées dans le graphique avec une seule décimale
            ax.scatter(max_x, max_y, color=color, zorder=10)
            ax.text(max_x + dx, max_y + dy, f'({max_x:.1f}, {max_y:.1f})', fontsize=14, fontweight='bold', ha=ha, va=va, color=color, bbox=dict(facecolor='white', edgecolor=color, boxstyle='round,pad=0.3'))
            plot_canvas.draw_idle()
        else:
            messagebox.showwarning("Warning", "Max point coordinates are out of range.")

    def choose_color():
        color = colorchooser.askcolor(color='red')[1]  # Choisissez la couleur par défaut ici
        if color:
            display_max_point(color)

    choose_color()
    
def add_data_point(ax, plot_canvas):
    # Demander à l'utilisateur de choisir une couleur
    color = colorchooser.askcolor()[1]
    if color:
        # Récupérer les données x et y du graphique
        x_data = ax.lines[0].get_xdata()
        y_data = ax.lines[0].get_ydata()

        # Ajouter les points avec la couleur spécifiée à l'aide de ax.scatter
        
        ax.scatter(x_data, y_data, color=color, zorder=10) 

        # Mettre à jour le canvas pour afficher les modifications
        plot_canvas.draw_idle()

def set_line_color(ax, plot_canvas):
    color = colorchooser.askcolor()[1]
    line_style = simpledialog.askstring("Line Style", "Enter line style ('-', '--', ':', '-.'):")
    if color and line_style:
        ax.lines[0].set_color(color)  # Change seulement la couleur de la première ligne (la courbe)
        ax.lines[0].set_linestyle(line_style)  
        plot_canvas.draw_idle()

def set_axes_color(ax, plot_canvas):
    color = colorchooser.askcolor()[1]
    if color:
        for spine in ax.spines.values():
            spine.set_color(color)
        
        # Changer la couleur des valeurs sur les axes
        ax.tick_params(axis='x', colors=color)
        ax.tick_params(axis='y', colors=color)

        plot_canvas.draw_idle()


def set_grid_color(ax, plot_canvas):
    color = colorchooser.askcolor()[1]
    if color:
        ax.grid(color=color)
        plot_canvas.draw_idle()

def set_label_color(ax, plot_canvas):
    color = colorchooser.askcolor()[1]
    if color:
        ax.xaxis.label.set_color(color)
        ax.yaxis.label.set_color(color)
        ax.title.set_color(color)
        plot_canvas.draw_idle()

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
    
def open_graph_settings_window(fig, plot_canvas):
    graph_settings_window = tk.Toplevel()
    graph_settings_window.title("Graph Settings")

    # Boutons pour les paramètres du graphique
    set_line_color_button = tk.Button(graph_settings_window, text="Set Line Color", command=lambda: set_line_color(fig.axes[0], plot_canvas))
    set_line_color_button.pack(side=tk.TOP, padx=5, pady=5)

    set_axes_color_button = tk.Button(graph_settings_window, text="Set Axis Color", command=lambda: set_axes_color(fig.axes[0], plot_canvas))
    set_axes_color_button.pack(side=tk.TOP, padx=5, pady=5)

    set_grid_color_button = tk.Button(graph_settings_window, text="Set grid", command=lambda: set_grid_color(fig.axes[0], plot_canvas))
    set_grid_color_button.pack(side=tk.TOP, padx=5, pady=5)

    set_label_color_button = tk.Button(graph_settings_window, text="Set Label Color", command=lambda: set_label_color(fig.axes[0], plot_canvas))
    set_label_color_button.pack(side=tk.TOP, padx=5, pady=5)

    set_background_color_button = tk.Button(graph_settings_window, text="Set Background Color", command=lambda: set_background_color(fig.axes[0], plot_canvas))
    set_background_color_button.pack(side=tk.TOP, padx=5, pady=5)

    # Bouton pour fermer la fenêtre des paramètres du graphique
    close_button = tk.Button(graph_settings_window, text="Close", command=graph_settings_window.destroy)
    close_button.pack(side=tk.BOTTOM, padx=5, pady=5)


def toggle_grid(ax, plot_canvas):
    # Set the grid lines visibility to True
    ax.grid(True)
    plot_canvas.draw_idle()

def set_scale(ax, plot_canvas):
    x_spacing = simpledialog.askfloat("X Spacing", "Enter spacing between x ticks:")
    y_spacing = simpledialog.askfloat("Y Spacing", "Enter spacing between y ticks:")
    if x_spacing and y_spacing:
        ax.xaxis.set_major_locator(plt.MultipleLocator(x_spacing))
        ax.yaxis.set_major_locator(plt.MultipleLocator(y_spacing))
        plt.tight_layout()
        plot_canvas.draw_idle()

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



def display_graph(fig):
    global graph_window
    if graph_window:
        graph_window.destroy()

    graph_window = tk.Toplevel()
    graph_window.title("Graph")
    graph_window.geometry("1000x600")


    plot_canvas = FigureCanvasTkAgg(fig, master=graph_window)
    plot_canvas.draw()
    plot_canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=1)

    toolbar = NavigationToolbar2Tk(plot_canvas, graph_window)
    toolbar.update()

    
    for ax in fig.get_axes():
        ax.set_autoscale_on(False)
        ax.set_xlim(auto=True)
        ax.set_ylim(auto=True)
        ax.lines[0].set_linewidth(1)  # Reset line width


    custom_button = tk.Button(graph_window, text="Custom Labels and Title", command=lambda: set_custom_labels_and_title(fig.axes[0], plot_canvas))
    custom_button.pack(side=tk.LEFT, padx=5)

    add_data_point_button = tk.Button(graph_window, text="Add Data Point", command=lambda: add_data_point(fig.axes[0], plot_canvas))
    add_data_point_button.pack(side=tk.LEFT, padx=5)

    set_scale_button = tk.Button(graph_window, text="Set Scale", command=lambda: set_scale(fig.axes[0], plot_canvas))
    set_scale_button.pack(side=tk.LEFT, padx=5)

    grid_button = tk.Button(graph_window, text="Toggle Grid", command=lambda: toggle_grid(fig.axes[0], plot_canvas))
    grid_button.pack(side=tk.LEFT, padx=5)
    
    display_max_button = tk.Button(graph_window, text="Display Max Point", command=lambda: display_max_point_and_coords(fig.axes[0], plot_canvas))
    display_max_button.pack(side=tk.LEFT, padx=5)

    graph_settings_button = tk.Button(graph_window, text="Graph Settings", command=lambda: open_graph_settings_window(fig, plot_canvas))
    graph_settings_button.pack(side=tk.LEFT, padx=5, pady=5)
    
    
    def on_closing_graph():
        if messagebox.askokcancel("Quit", "Are you sure you want to quit?"):
            graph_window.withdraw()  # Hide the window
            graph_window.quit()  # Quit the main loop

            
    graph_window.protocol("WM_DELETE_WINDOW", on_closing_graph)
    graph_window.mainloop()

