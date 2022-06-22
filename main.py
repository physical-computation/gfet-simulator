# Main file for GFET Lab
import GUI as gui
import tkinter as tk

def main():
    root = tk.Tk()
    app = gui.GUI(root)
    root.mainloop()

if __name__ == '__main__':
    main()
