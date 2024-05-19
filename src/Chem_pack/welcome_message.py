import tkinter as tk
window = None

def welcome_message():
    global entry_input
    welcome_window = tk.Toplevel(window)
    welcome_window.title("Welcome Message")

    welcome_text = (
        "✨ Welcome to Our Project! ✨\n\n"
        "We are excited to share our project with you.\n\n"
        "Discover all the details on our GitHub repository:\n"
        "🔗 [Project ppchem Tools Kit] (https://github.com/sgrunber/Project-ppchem-tools-kit)\n\n"
        "Happy exploring! 🚀\n\n"
    )

    welcome_textbox = tk.Text(welcome_window, wrap="word", font=("Times New Roman", 20), fg="#BBE1FA", bg = "#1B262C", height=10, width=50)
    welcome_textbox.insert("1.0", welcome_text)
    welcome_textbox.config(state="disabled")
    welcome_textbox.tag_configure("center", justify="center")
    welcome_textbox.tag_add("center", "1.0", "end")
    welcome_textbox.grid(row=0, column=0, padx=0, pady=0)

