import numpy as np
from tkinter import messagebox, colorchooser


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
            dx = 0.03 * (x_lim[1] - x_lim[0])
            dy = 0.03 * (y_lim[1] - y_lim[0])
            
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