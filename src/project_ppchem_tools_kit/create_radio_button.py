import tkinter as tk
from tkinter import Canvas
window = None
selected_radio=None
from project_ppchem_tools_kit.on_radio_select import on_radio_select

canvas = Canvas(window, bg="#1B262C", height=800, width=1000, bd=0, highlightthickness=0, relief="ridge")



def create_radio_button(x, y, text, value):
    """Creates a radio button with the specified text and value at the specified position (x, y)

    Args:
        x (int): The x-coordinate for placing the radio button.
        y (int): The y-coordinate for placing the radio button.
        text (str): The label text for the radio button.
        value (str): The value associated with the radio button.
    """
    radio_button = tk.Radiobutton(canvas, text=text, variable=selected_radio, value=value,
                                  command=lambda: on_radio_select(value),
                                  font=("Times New Roman", 20), bg="#1B262C")
    canvas.create_window(x, y, anchor="nw", window=radio_button)
