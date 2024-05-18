import tkinter as tk
window = None

def welcome_message():
    global entry_input
    welcome_window = tk.Toplevel(window)
    welcome_window.title("Welcome Message")

    welcome_text = (
        "Welcome!\n\n"
        "Here is our project: https://github.com/sgrunber/Project-ppchem-tools-kit\n\n"
        "Enjoy ;)"
    )
    welcome_textbox = tk.Text(welcome_window, wrap="word", font=("Times New Roman", 25), fg="#BBE1FA", bg = "#1B262C", height=7, width=45)
    welcome_textbox.insert("1.0", welcome_text)
    welcome_textbox.config(state="disabled")
    welcome_textbox.grid(row=0, column=0, padx=0, pady=0)

