from project_ppchem_tools_kit.process_input import process_input
window = None

def bind_enter(event):
    try:
        if event.keysym == "Return":
            process_input()
    except Exception as e:
        print("Error while binding Enter key:", e)

# Attachez le gestionnaire d'événements à la fenêtre principale
if window:
    window.bind("<KeyPress>", bind_enter)