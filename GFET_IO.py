# Module to handle all I/O Functionality of GFET Simulator

import tkinter as tk
from tkinter import filedialog

import numpy as np

class GFET_IO:

    def __init__(self, *args, **kwargs):
        self.data = {}

    # Load Values (dielectric materials etc. from txt file
    def loadDielectrics():
        dielectrics = np.loadtxt('Dielectrics.txt', dtype=np.dtype('O'), delimiter=',', skiprows=1)
        return dielectrics

    def exportData(self, data):
        f = filedialog.asksaveasfile(mode='w', defaultextension=".csv")
        if f is None:
            return
    
        # Group voltages and currents into pairs. For general data exporting,
        # maybe some unified structure needed...
        self.data.update(data)
        
        dataToSave = ''.join([str(i)+","+str(j)+"\n" for i,j in zip(self.data["Vgs"],self.data["Ids"])])
    
        f.write(dataToSave)
        f.close()

    def loadSweep(self):
        f = filedialog.askopenfile(mode='r', filetypes=[('CSV Files', '*.csv')])

        if f is not None:
            content = f.read().split('\n')
            
            self.data["V"] = [float(item) for item in content]
            extSweep = True
            datapoints = len(self.data["V"])
