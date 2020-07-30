# Module to handle all I/O Functionality of GFET Simulator

# For Physical Constants
from scipy import constants as consts

import tkinter as tk
from tkinter import filedialog

import numpy as np
import itertools
import csv

class GFET_IO:

    def __init__(self, *args, **kwargs):
        self.transData = {}
        self.ivData = {}
        self.extTransSweep = False
        self.extIVSweep = False

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

    def exportFreq(self, data):
        filename = filedialog.asksaveasfilename()#mode='w', defaultextension=".csv")
        if filename is None:
            return
    
        self.transData.update(data["TransChars"])
        
        dataPairs = []
        

        for index, Id in enumerate(self.transData["Ids"]): # for each entry in Ids
            column = []
            for datapoint in range(len(self.transData["Ids"][index])): #for each datapoint in that index
                column.append(str(self.transData["Vtg"][datapoint]) + ','
                              + str(self.transData["Ids"][index][datapoint]) + ','
                              + str(data["fT"][index][datapoint]))
            dataPairs.append(column)

        rows = list(zip(*itertools.chain(dataPairs)))

        headerRow = []
        for dp in self.transData["Vds"]:
            row = "Vds:" + ',' + str(dp)
            headerRow.append(row)

        titlerow = []
        for dp in self.transData["Vds"]:
            row = "Vtg (V):" + "," + "Ids (A):" + "," + "fT (Hz)"
            titlerow.append(row)
        
        with open(filename + '.csv', 'w') as f:
            writer = csv.writer(f, quoting=csv.QUOTE_NONE, escapechar=" ")
            writer.writerow(headerRow)
            writer.writerow("") # Blank spacer row, for formatting
            writer.writerow(titlerow)
            for row in rows:
                writer.writerow(row)
            f.close()

    def loadSweep(self, sweepType):
        try:
            f = filedialog.askopenfilename(filetypes=[('CSV Files', '*.csv')])
        except FileNotFoundError:
            return

        if f == '':
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

            if sweepType == "Gate":
                for column in range(columns):
                    Vds.append(float(row1[column]))
            elif sweepType == "Drain":            
                for column in range(columns):
                    Vtg.append(float(row1[column]))

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

            # Vbg, if applicable. Stub for now
            Vbg = 0.0


            # Add first value to sweep. scientifically weird to
            # change both bias and gate sweep, so assumes only one
            # set of gate sweep values...
            if sweepType == "Gate":
                Vtg.append(float(row1[0]))
            elif sweepType == "Drain":
                Vds.append(float(row1[0]))
                
            # Now get the sweeps. Should be third row onwards
            if sweepType == "Gate":
                for row in reader2:
                    Vtg.append(float(row[0]))
            elif sweepType == "Drain":                    
                for row in reader2:
                    Vds.append(float(row[0]))
                    
            # Finally, update sweep data
            if sweepType == "Gate":
                self.extTransSweep = True
                self.transData.update({"Vds": Vds,
                                       "Vtg": Vtg,
                                       "Vbg": Vbg})
            elif sweepType == "Drain":
                self.extIVSweep = True
                self.transData.update({"Vds": Vds,
                                       "Vtg": Vtg,
                                       "Vbg": Vbg})

    def expTemp(self, biasVoltage, sweepVoltage):
        filename = filedialog.asksaveasfilename() #Maybe Give a default name, but choose location
        if filename is None:
            return
        
        # First few rows as an example:        
        with open(filename + '.csv', 'w') as f:
            writer = csv.writer(f, quoting=csv.QUOTE_NONE, escapechar=" ")
            writer.writerow([biasVoltage])
            writer.writerow(["x,x,..."])
            writer.writerow("")
            writer.writerow("") # Blank spacer row, for formatting
            writer.writerow([sweepVoltage])
            writer.writerow("x")
            writer.writerow("x")
            writer.writerow(["..."])
            f.close()
            
    def exportSPICEModel(self, model, params, eps):
        filename = filedialog.asksaveasfilename()
        if filename is None:
            return

        # Parse model parameters from variables
        parameters = []
        
        parameters.append("+L= " + str(float(params[3])*10**(-6)))
        parameters.append("+W=" + str(float(params[2])*10**(-6)))
        parameters.append("+Tox=" + str(float(params[0])*10**(-9)))
        parameters.append("+er= " + str(eps.get().split("(")[1].replace(")","")))
        parameters.append("+mu= " + str(float(params[4])*10**(-4))) # cm2/Vs to m2/Vs
        parameters.append("+omega= " + str(float(params[5])/consts.hbar))
        parameters.append("+Nf= " + str(float(params[6])))

        # Write constants
        parameters.append("+e= " + str(consts.elementary_charge))
        parameters.append("+pi= " + str(consts.pi))
        parameters.append("+e0= " + str(consts.epsilon_0))

        # Probably will want to make the models a bit more dynamic
        # in terms of loads, rather than hard-coding
        if model == "Rodriguez":
            headerLine = "*                      G  Dc Do S"
            subcktLine = ".subckt GFET_Rodriguez n1 n2 n3 n4"
            paramsLine = ".params"
            modelDef = "B1 n1 n3 I = {abs((mu*W*Ctg*((V(n1)+Vth)-V(n2,n4)/2))/((L/V(n2,n4))+(mu/omega)*sqrt((pi*Ctg)/e)*sqrt((V(n1)+Vth)-V(n2,n4)/2)))}"
            endsLine = ".ends GFET_Rodriguez"
                          
            with open(filename + '.lib', 'w') as f:
                f.write(headerLine + "\n")
                f.write(subcktLine + "\n")
                f.write(paramsLine + "\n")

                # write each parameter
                for parameter in parameters:
                    f.write(parameter + "\n")

                # Write calculations
                f.write("+eox = {er*e0}\n")
                f.write("+Ctg = {eox/Tox}\n")
                f.write("+Vth = {e*Nf/Ctg}\n")
                
                # Write model definition
                f.write(modelDef + "\n")
                f.write(endsLine + "\n")
                f.close()
            print("model successfully output")
