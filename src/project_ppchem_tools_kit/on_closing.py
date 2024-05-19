from tkinter import messagebox
window=None

def on_closing():
    global window
    if messagebox.askokcancel("Quit", "Are you sure you want to quit?"):
        window.destroy()
window.protocol("WM_DELETE_WINDOW", on_closing)
