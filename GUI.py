# Class to setup and manage GUI components for GFET simulator program.

import GFET_IO as gio
import models as gfet

import tkinter as tk
from tkinter import ttk, filedialog

import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import (
    FigureCanvasTkAgg, NavigationToolbar2Tk)
from matplotlib.backend_bases import key_press_handler
from matplotlib.figure import Figure
from matplotlib.widgets import Cursor

import numpy as np

#from mpldatacursor import datacursor

# Device Parameters
tox1 = 50 # TG oxide layer thickness, metres
W = 40 # Channel Width, metres
L = 120 # Channel Length, metres
mu = 7000 # Effective carrier mobility
N = 0.5*10**(16) # carrier density
Ep = 0.56*10**(-19) # Surface phonon energy of the substrate
Vgt0 = -0.8229 # TG Voltage at Dirac pt, volts. (Avg from measurements)

datapoints = 100

fields1 = '-Vgs', '+Vgs', 'Vds', 'Datapoints'
params1 = -10, 10, 0.2, datapoints

fields2 = 'Dielectric Thickness (nm)', 'Channel Width (um)', 'Channel Length (um)', 'Mobility', 'Phonon Energy', 'Carrier Density', 'Gate Threshold\nVoltage'
params2 = tox1, W, L, mu, Ep, N, Vgt0

default_resolution = [800, 600]

class GUI:

    def __init__(self, master):
        self.root = master
        self.root.geometry(str(default_resolution[0]) + 'x' + str(default_resolution[1]))
        self.root.title('GFET Simulation')

        self.ents1 = ''
        self.ents2 = ''
        self.model = None
        self.ax = None
        self.canvas = None

        self.modelCombo = None
        self.dielecCombo = None
        self.data = {}
        
        left = tk.Frame(self.root, width=int(default_resolution[0]/2), height=default_resolution[1])
        right = tk.Frame(self.root, width=int(default_resolution[0]/2), height=default_resolution[1])

        # layout all of the main containers
        self.root.grid_rowconfigure(1, weight=1)
        self.root.grid_columnconfigure(0, weight=1)

        left.pack(side="left")
        right.pack(side="right")

        # Left Panel Tabs (settings etc)
        tab_control = ttk.Notebook(left)
        tab_control.pack()

        self.tab1 = ttk.Frame(tab_control)
        tab_control.add(self.tab1,text= 'Sweep Parameters')
                        
        self.tab2 = ttk.Frame(tab_control)
        tab_control.add(self.tab2, text= 'Device Parameters')

        self.setupParamsTab(self.tab2, self.root)
        self.setupSweepTab(self.tab1, self.root)

        # Right Panel Tabs (plots etc)
        tab_control2 = ttk.Notebook(right)
        tab_control2.pack()

        self.tab1 = ttk.Frame(tab_control2)
        tab_control2.add(self.tab1,text= 'Transfer Characteristic')
                        
        self.tab2 = ttk.Frame(tab_control2)
        tab_control2.add(self.tab2, text= 'I-V Characteristic')
        self.setupAxes(right)

    def fetch(self, entries):
        data = []
        for entry in entries:
            text  = entry[1].get()
            data.append(text)
        return data

    def makeform(self, root, fields, params):
        entries = []
        for index,field in enumerate(fields):
            row = tk.Frame(root)
            lab = tk.Label(row, width=15, text=field, anchor='w')

            ent = tk.Entry(row, bd=1)
            ent.insert(0, str(params[index]))
            
            row.pack(side=tk.TOP, fill=tk.X, padx=5, pady=5)
            lab.pack(side=tk.LEFT)
            ent.pack(side=tk.RIGHT, expand=tk.YES, fill=tk.X)
            entries.append((field, ent))
        return entries

    def generateSweep(self, ent):

        ent1 = self.fetch(ent)

        start_voltage = float(ent1[0])
        end_voltage = float(ent1[1])
        Vds = float(ent1[2])
        dps = int(ent1[3])
        
        Vgt = list(np.linspace(start_voltage, end_voltage, dps))

        return {"Vgs": Vgt,
                "Vds": Vds,
                "dps": dps}

    def loadModel(self, name, ents, ents2):
        self.model = name

        # Todo: figure out loading sweeps vs getting from the GUI
        self.data.update(self.generateSweep(ents))

        params = self.fetch(ents2)
        eps = self.dielecCombo

        if self.model == 'Rodriguez':
            GFET = gfet.RodriguezGFET(params, self.data["Vds"], self.data["Vgs"], self.data["dps"], eps)
            self.data.update(GFET.calculateIds())
        elif self.model == 'Thiele':
            GFET = gfet.ThieleGFET(params, self.data["Vds"], self.data["Vgs"], self.data["dps"], eps)
            self.data.update(GFET.calculateIds())
        else:
            print('ERROR: No models available.')

        if self.model != None:
            self.plotCurve() # maybe all curves?

    def plotCurve(self):
        self.ax.clear()
        self.ax.set_title('Transfer Characteristic')
        self.ax.set_ylabel(r'$I_{DS}$ (A)')
        self.ax.set_xlabel(r'$V_{GS}$ (V)')
        lines = self.ax.plot(self.data["Vgs"], self.data["Ids"])
        self.ax.set_aspect(1./self.ax.get_data_ratio())
        self.ax.ticklabel_format(axis='y',style='sci', scilimits=(0,0))    
    #    datacursor(lines, display='single')
        self.canvas.draw()       
    
    def setupSweepTab(self, tab1, root):
        
        self.ents = self.makeform(tab1, fields1, params1)
        root.bind('<Return>', (lambda event, e=self.ents: self.fetch(e)))

        modelLabel = tk.Label(tab1, text='GFET Model')
        self.modelCombo = ttk.Combobox(tab1, state='readonly')

        self.modelCombo['values'] = 'Rodriguez', 'Thiele'
        self.modelCombo.current(0) # First value in list is default
        
        root.bind("<<ComboboxSelected>>")

        modelLabel.grid(row=1, column=0)
        self.modelCombo.grid(row=1,column=1)
        modelLabel.pack()
        self.modelCombo.pack()

        io = gio.GFET_IO()

        b1 = tk.Button(tab1, text='Simulate',
                  command=(lambda e=self.ents, e2=self.ents2: self.loadModel(self.modelCombo.get(), e, e2)))
        b1.pack(side=tk.LEFT, padx=5, pady=5)
        b2 = tk.Button(tab1, text='Load Sweep', command=io.loadSweep)
        b2.pack(side=tk.LEFT, padx=5, pady=5)
        b3 = tk.Button(tab1, text='Export CSV', command=(lambda : io.exportData(self.data)))
        b3.pack(side=tk.LEFT, padx=5, pady=5)
        
    def setupParamsTab(self, tab2, root):
        
        self.ents2 = self.makeform(tab2, fields2, params2)
        root.bind('<Return>', (lambda event, e=self.ents2: self.fetch(e)))

        dielecLabel = tk.Label(tab2, text='Relative Permittivity')
        self.dielecCombo = ttk.Combobox(tab2, state='readonly')

        dielecs = gio.GFET_IO.loadDielectrics()
        
        vals = []
        i = 0
        for i in range(len(dielecs)):
            vals.append(dielecs[i][0] + " (" + str(dielecs[i][1]) + ")")

        self.dielecCombo['values'] = vals
        self.dielecCombo.current(0) # First value in list is default
        
        root.bind("<<ComboboxSelected>>")
        dielecLabel.grid(row=1, column=0)
        self.dielecCombo.grid(row=1,column=1)
        dielecLabel.pack()
        self.dielecCombo.pack()

    def setupAxes(self, frame):
        fig = plt.Figure(figsize=(5,4), dpi=100)
        self.ax = fig.add_subplot(111)
        self.ax.set_title('Transfer Characteristic')
        self.ax.set_ylabel(r'$I_{DS}$ (A)')
        self.ax.set_xlabel(r'$V_{GS}$ (V)')
        self.ax.set_aspect(1./self.ax.get_data_ratio())
        self.ax.ticklabel_format(axis='both',style='sci', scilimits=(0,0))

        self.canvas = FigureCanvasTkAgg(fig, master=frame)
        self.canvas.get_tk_widget().pack()
        self.canvas.draw()
