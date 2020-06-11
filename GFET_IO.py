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
        f = filedialog.askopenfilename(filetypes=[('CSV Files', '*.csv')])

        if f is None:
            return
        
        # Get the bias voltage(s) first
        with open(f, newline='') as csvfile:
            reader1, reader2 = itertools.tee(csv.reader(csvfile, delimiter=','))

            columns = len(next(reader1))
            del reader1

            # Set up lists for each variable
            Vds = []
            Vtg = []

            # Get biases first. Always second row
            next(reader2)
            row1 = next(reader2)
            
            for column in range(columns):
                Vds.append(float(row1[column]))

            sweep = False
            
            while sweep == False:
                row1 = next(reader2) # iterate row in csv
            
                # if the first number can be converted to a float,
                # i.e. is a number, not a heading, break loop
                try: 
                    float(row1[0])
                    sweep = True
                except (ValueError, IndexError) as e:
                    sweep = False

            # Add first value to sweep. scientifically weird to
            # change both bias and gate sweep, so assumes only one
            # set of gate sweep values...
            Vtg.append(float(row1[0]))
                
            # Now get the sweeps. Should be third row onwards
            for row in reader2:
                Vtg.append(float(row[0]))

            print("\n Values Extracted: ")
            print(Vds)
            print(Vtg)
        
            # Finally, update sweep data
            self.transData.update({"Vds": Vds})
            self.transData.update({"Vtg": Vtg})
            self.extSweep = True
                
