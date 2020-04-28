# Main file for GFET Simulator

import GUI as gui

import tkinter as tk

##class Main(tk.Frame):
##
##    def __init__(self, parent, *args, **kwargs):
##        tk.Frame.__init__(self, parent, *args, **kwargs)
##        self.parent = parent

def main():
    root = tk.Tk()
    app = gui.GUI(root)
    root.mainloop()

if __name__ == '__main__':
    main()
