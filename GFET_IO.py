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
            
    def exportSPICEModel(self, model, params, eps1, eps2):
        if model == "Rodriguez":
            filename = filedialog.asksaveasfilename()
            if filename is None:
                return

            # Parse model parameters from variables
            parameters = []
            
            parameters.append("+L= " + str(float(params[3])*10**(-6)))
            parameters.append("+W=" + str(float(params[2])*10**(-6)))
            parameters.append("+Tox=" + str(float(params[0])*10**(-9)))
            parameters.append("+er= " + str(eps1.get().split("(")[1].replace(")","")))
            parameters.append("+mu= " + str(float(params[4]))) # cm2/Vs to m2/Vs
            parameters.append("+omega= " + str((float(params[5])*10**(-19))/(2*consts.pi*consts.h)))
            parameters.append("+Nf= " + str(float(params[6])))

            # Probably will want to make the models a bit more dynamic
            # in terms of loads, rather than hard-coding
            
            headerLine = "*                      G  D  S"
            subcktLine = ".subckt GFET_Rodriguez n1 n2 n3"
            paramsLine = ".params"
            modelDef = "BI n1 n2 I = abs((mu*W*Ctg*((V(n1)+Vth-V(n2)/2)))/((L/V(n2))+(mu/omega)*sqrt((pi*Ctg)/echarge)*sqrt(abs(V(n1)+Vth-V(n2)/2))))"
            endsLine = ".ends GFET_Rodriguez"
                          
            with open(filename + '.lib', 'w') as f:
                f.write(headerLine + "\n")
                f.write(subcktLine + "\n")
                f.write(paramsLine + "\n")

                # write each parameter
                for parameter in parameters:
                    f.write(parameter + "\n")

                # Write constants
                f.write("+e0 = 8.854e-12\n")
                f.write("+pi = 3.141\n")
                
                # Write calculations
                f.write("+eox = {er*e0}\n")
                f.write("+Ctg = {eox/Tox}\n")
                f.write("+Vth = {(echarge*Nf)/Ctg}\n")
                f.write("+omega={hw/(2*pi*planck)}\n")
                
                # Write model definition
                f.write(modelDef + "\n")
                f.write(endsLine + "\n")
                f.close()
            tk.messagebox.showinfo("Note:", "Successfully exported SPICE model!")
        elif model == "Jimenez":
            filename = filedialog.asksaveasfilename()
            if filename is None:
                return

            # Parse model parameters from variables
            parameters = []
            
            parameters.append("+L= " + str(float(params[3])*10**(-6)))
            parameters.append("+W=" + str(float(params[2])*10**(-6)))
            parameters.append("+Tox1=" + str(float(params[0])*10**(-9)))
            parameters.append("+Tox2=" + str(float(params[1])*10**(-9)))
            parameters.append("+er1= " + str(eps1.get().split("(")[1].replace(")","")))
            parameters.append("+er2= " + str(eps1.get().split("(")[1].replace(")","")))
            parameters.append("+mu= " + str(float(params[4]))) # cm2/Vs to m2/Vs
            parameters.append("+hw= {" + str((float(params[5]))) + "*echarge}")
            parameters.append("+Nf= " + str(float(params[6])))

            # Probably will want to make the models a bit more dynamic
            # in terms of loads, rather than hard-coding
            
            headerLine = "*                    G  D  S"
            subcktLine = ".subckt GFET_Jimenez n1 n2 n3"
            paramsLine = ".params"
            modelDef = "BI n1 n2 I = ((mu*k)/2)*(W/Leff(V(n2)))*(g1(V(n1),V(n2))-g2(V(n1)))"
            endsLine = ".ends GFET_Jimenez"
                          
            with open(filename + '.lib', 'w') as f:
                f.write(headerLine + "\n")
                f.write(subcktLine + "\n")
                f.write(paramsLine + "\n")

                # write each parameter
                for parameter in parameters:
                    f.write(parameter + "\n")

                f.write("+vF=1e6\n")
                f.write("+e0=8.854e-12\n")
                f.write("+Vbg=0.0\n")
                f.write("+Vg0=0.85\n")
                f.write("+Vb0=0\n")
                
                # Write calculations
                f.write("+hbar={planck/(2*pi)}\n")
                f.write("+omega={hw/hbar}\n")
                f.write("+eox1={er1*e0}\n")
                f.write("+eox2={er2*e0}\n")
                f.write("+Ct={eox1/Tox1}\n")
                f.write("+Cb={eox2/Tox2}\n")
                f.write("+k={((2*(echarge**2))/pi)*(echarge/((hbar*vF)**2))}\n")

                # Write functions
                f.write(".func Leff(vds) {L+mu*(abs(vds)/vF)}\n")
                f.write(".func Vcd(v,vds) {if((((v-Vg0-vds)*Ct+(Vbg-Vb0-vds)*Cb) > 0), ((-(Ct+Cb)+sqrt(((Ct+Cb)**2)+2*k*((v-Vg0-vds)*Ct+(Vbg-Vb0-vds)*Cb)))/(k)), ((-(Ct+Cb)+sqrt(((Ct+Cb)**2)-2*k*((v-Vg0-vds)*Ct+(Vbg-Vb0-vds)*Cb)))/(-k)))}\n")
                f.write(".func Vcs(v) {if((((v-Vg0)*Ct+(Vbg-Vb0)*Cb) > 0), ((-(Ct+Cb)+sqrt(((Ct+Cb)**2)+2*k*((v-Vg0)*Ct+(Vbg-Vb0)*Cb)))/(k)), ((-(Ct+Cb)+sqrt(((Ct+Cb)**2)-2*k*((v-Vg0)*Ct+(Vbg-Vb0)*Cb)))/(-k)))}\n")
                f.write(".func g1(v,vds) {(-(Vcd(v,vds)**3)/3)-sgn(Vcd(v,vds))*((k*(Vcd(v,vds)**4))/(4*(Ct+Cb)))}\n")
                f.write(".func g2(v) {(-(Vcs(v)**3)/3)-sgn(Vcs(v))*((k*(Vcs(v)**4))/(4*(Ct+Cb)))}\n")
                        
                # Write model definition
                f.write(modelDef + "\n")
                f.write(endsLine + "\n")
                f.close()
            tk.messagebox.showinfo("Note:", "Successfully exported SPICE model!")
        else:
            tk.messagebox.showinfo("Warning", "Model could not be created: not a SPICE-compatible model.")
