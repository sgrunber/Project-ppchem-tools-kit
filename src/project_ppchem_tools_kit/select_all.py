window = None

def select_all(event):
    if window and window.winfo_exists():
        event.widget.tag_add("sel", "1.0", "end")
    return "break"