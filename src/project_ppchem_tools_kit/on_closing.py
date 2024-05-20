from tkinter import messagebox
window=None

def on_closing():
    """
    Handle the closing event for the graph window.

    This function is triggered when the user attempts to close the graph window. 
    It prompts the user with a confirmation dialog to ensure they want to quit. 
    If the user confirms (clicks 'OK'), the graph window is hidden and the 
    main event loop is terminated.
    """
    global window
    if messagebox.askokcancel("Quit", "Are you sure you want to quit?"):
        window.destroy()
window.protocol("WM_DELETE_WINDOW", on_closing)
