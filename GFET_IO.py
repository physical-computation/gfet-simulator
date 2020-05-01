# Module to handle all I/O Functionality of GFET Simulator

import tkinter as tk
from tkinter import filedialog

import numpy as np
import itertools
import csv

class GFET_IO:

    def __init__(self, *args, **kwargs):
        self.data = {}

    # Load Values (dielectric materials etc. from txt file
    def loadDielectrics():
        dielectrics = np.loadtxt('Dielectrics.txt', dtype=np.dtype('O'), delimiter=',', skiprows=1)
        return dielectrics

    def exportTransferChars(self, data):
        filename = filedialog.asksaveasfilename()#mode='w', defaultextension=".csv")
        if filename is None:
            return
    
        # Group voltages and currents into pairs. For general data exporting,
        # maybe some unified structure needed...
        self.data.update(data)
        
        dataPairs = []

        for index, Id in enumerate(self.data["Ids"]): # for each entry in Ids
            column = []
            for datapoint in range(len(self.data["Ids"][index])): #for each datapoint in that index
                column.append(str(self.data["Vgs"][datapoint]) + ',' + str(self.data["Ids"][index][datapoint]))
            dataPairs.append(column)

        rows = list(zip(*itertools.chain(dataPairs)))

        headerRow = []
        for dp in self.data["Vds"]:
            row = "Vds:" + ',' + str(dp)
            headerRow.append(row)

        titlerow = []
        for dp in self.data["Vds"]:
            row = "Vgs (V):" + "," + "Ids (A):"
            titlerow.append(row)
        
        with open(filename + '.csv', 'w') as f:
            writer = csv.writer(f, quoting=csv.QUOTE_NONE, escapechar=" ")
            writer.writerow(headerRow)
            writer.writerow("") # Blank spacer row, for formatting
            writer.writerow(titlerow)
            for row in rows:
                writer.writerow(row)
            f.close()

    def exportIVChars(self, data):
        filename = filedialog.asksaveasfilename()#mode='w', defaultextension=".csv")
        if filename is None:
            return
    
        # Group voltages and currents into pairs. For general data exporting,
        # maybe some unified structure needed...
        self.data.update(data)

        print("\nIV Chars Vgs:")
        print(self.data["Vgs"])

        print("\n\nIV Chars Ids:")
        print(self.data["Ids"])

        print("\n\nIV Chars Vds:")
        print(self.data["Vds"])
        
##        dataPairs = []
##
##        for index, Id in enumerate(self.data["Ids"]): # for each entry in Ids
##            column = []
##            for datapoint in range(len(self.data["Ids"][index])): #for each datapoint in that index
##                column.append(str(self.data["Vds"][datapoint]) + ',' + str(self.data["Ids"][index][datapoint]))
##            dataPairs.append(column)
##
##        rows = list(zip(*itertools.chain(dataPairs)))
##        print(dataPairs)
##        headerRow = []
##        for dp in self.data["Vgs"]:
##            row = "Vgs:" + ',' + str(dp)
##            headerRow.append(row)
##
##        titlerow = []
##        for dp in self.data["Vgs"]:
##            row = "Vds (V):" + "," + "Ids (A):"
##            titlerow.append(row)
##        
##        with open(filename + '.csv', 'w') as f:
##            writer = csv.writer(f, quoting=csv.QUOTE_NONE, escapechar=" ")
##            writer.writerow(headerRow)
##            writer.writerow("") # Blank spacer row, for formatting
##            writer.writerow(titlerow)
##            for row in rows:
##                writer.writerow(row)
##            f.close()

    def loadSweep(self):
        f = filedialog.askopenfile(mode='r', filetypes=[('CSV Files', '*.csv')])

        if f is not None:
            content = f.read().split('\n')
            
            self.data["V"] = [float(item) for item in content]
            extSweep = True
            datapoints = len(self.data["V"])
