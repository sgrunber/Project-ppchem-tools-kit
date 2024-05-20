window = None

def select_all(event):
    """Selects all text in a widget when an event is triggered.

    Args:
        event: The event that triggered the function.

    Returns:
        str: "break" to prevent further propagation of the event.
    """
    if window and window.winfo_exists():
        event.widget.tag_add("sel", "1.0", "end")
    return "break"