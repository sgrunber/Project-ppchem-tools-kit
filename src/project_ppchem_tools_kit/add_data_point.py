from tkinter import colorchooser

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