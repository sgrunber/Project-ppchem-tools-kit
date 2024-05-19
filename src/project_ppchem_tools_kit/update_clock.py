window = None
import tkinter as tk
from datetime import datetime

def update_clock():
    current_time = datetime.now().strftime('%H:%M:%S')
    current_date = datetime.now().strftime('%Y-%m-%d')
    clock_label.config(text=current_time)
    date_label.config(text=current_date)
    window.after(1000, update_clock)

date_label = tk.Label(window, text="", font=("Times New Roman", 20), bg="#1B262C", fg="#FFFFFF")
date_label.place(x=890, y=10)

clock_label = tk.Label(window, text="", font=("Times New Roman", 20), bg="#1B262C", fg="#FFFFFF")
clock_label.place(x=890, y=40)