# Module to handle all I/O Functionality of GFET Simulator

import tkinter as tk
from tkinter import filedialog

import numpy as np
import itertools
import csv

class GFET_IO:

    def __init__(self, *args, **kwargs):
        self.transData = {}
        self.ivData = {}
        self.extSweep = False

    # Load Values (dielectric materials etc. from txt file
    def loadDielectrics(self):
        dielectrics = np.loadtxt('Dielectrics.txt', dtype=np.dtype('O'), delimiter=',', skiprows=1)
        return dielectrics

    def exportTransferChars(self, data):
        filename = filedialog.asksaveasfilename()#mode='w', defaultextension=".csv")
        if filename is None:
            return
    
        self.transData.update(data["TransChars"])
        
        dataPairs = []

        for index, Id in enumerate(self.transData["Ids"]): # for each entry in Ids
            column = []
            for datapoint in range(len(self.transData["Ids"][index])): #for each datapoint in that index
                column.append(str(self.transData["Vtg"][datapoint]) + ','
                              + str(self.transData["Ids"][index][datapoint]))
            dataPairs.append(column)

        rows = list(zip(*itertools.chain(dataPairs)))

        headerRow = []
        for dp in self.transData["Vds"]:
            row = "Vds:" + ',' + str(dp)
            headerRow.append(row)

        titlerow = []
        for dp in self.transData["Vds"]:
            row = "Vtg (V):" + "," + "Ids (A):"
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
        
        self.ivData.update(data["IVChars"])
        
        dataPairs = []

        for index, Id in enumerate(self.ivData["Ids"]): # for each entry in Ids
            column = []
            for datapoint in range(len(self.ivData["Ids"][index])): #for each datapoint in the index
                column.append(str(self.ivData["Vds"][datapoint]) + ','
                              + str(self.ivData["Ids"][index][datapoint]))
            dataPairs.append(column)

        rows = list(zip(*itertools.chain(dataPairs)))
        
        headerRow = []
        for dp in self.ivData["Vtg"]:
            row = "Vtg:" + ',' + str(dp)
            headerRow.append(row)

        titlerow = []
        for dp in self.ivData["Vtg"]:
            row = "Vds (V):" + "," + "Ids (A):"
            titlerow.append(row)
        
        with open(filename + '.csv', 'w') as f:
            writer = csv.writer(f, quoting=csv.QUOTE_NONE, escapechar=" ")
            writer.writerow(headerRow)
            writer.writerow("") # Blank spacer row, for formatting
            writer.writerow(titlerow)
            for row in rows:
                writer.writerow(row)
            f.close()

    def loadSweep(self):
        f = filedialog.askopenfile(mode='r', filetypes=[('CSV Files', '*.csv')])

        if f is not None:
            content = f.read().split('\n')
            
            data = [float(item) for item in content]

            self.transData.update({"Vtg": data})
            self.extSweep = True
            
