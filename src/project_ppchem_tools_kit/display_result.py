import tkinter as tk
window = None

def display_result(result_text):
    """Displays the result of an operation in a new window

    Args:
        result_text (str): The text to be displayed in the result window
    """
    result_window = tk.Toplevel(window)
    result_window.title("Result")


    result_textbox = tk.Text(result_window, wrap="word", font=("Times New Roman", 25), fg="#BBE1FA", bg = "#1B262C", height=10, width=30)
    result_textbox.insert("1.0", result_text)
    result_textbox.config(state="disabled")
    result_textbox.grid(row=0, column=0, padx=0, pady=0)
