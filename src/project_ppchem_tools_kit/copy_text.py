window = None 

def copy_text(event):
    """Copies selected text to the clipboard. Bound to an event.

    Args:
        event: The event that triggered the function.

    Returns:
        None: Does not return any values.
    """
    try:
        if window and window.winfo_exists():
            if event and event.widget and event.widget.winfo_exists():
                event.widget.event_generate("<<Copy>>")
    except Exception as e:
        print("Error while copying text:", e)
    return "break"