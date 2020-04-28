# Module to handle all I/O Functionality of GFET Simulator

import tkinter as tk
from tkinter import filedialog

import numpy as np

class GFET_IO:

    def __init__(self):
        self.V = []
        self.I = []

    # Load Values (dielectric materials etc. from txt file
    def loadDielectrics():
        dielectrics = np.loadtxt('Dielectrics.txt', dtype=np.dtype('O'), delimiter=',', skiprows=1)
        return dielectrics

    def exportData(self, V, I):
        f = filedialog.asksaveasfile(mode='w', defaultextension=".csv")
        if f is None:
            return

        # Group voltages and currents into pairs
        dataToSave = ''.join([str(i)+","+str(j)+"\n" for i,j in zip(self.V,self.I)])
    
        f.write(dataToSave)
        f.close()

    def loadSweep(self):
        f = filedialog.askopenfile(mode='r', filetypes=[('CSV Files', '*.csv')])

        if f is not None:
            content = f.read().split('\n')
            
            self.V = [float(item) for item in content]
            print(V)
            extSweep = True
            datapoints = len(self.V)
