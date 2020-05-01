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
N = 5e17 # carrier density
Ep = 5.6e-20 # Surface phonon energy of the substrate
Vgt0 = -0.8229 # TG Voltage at Dirac pt, volts. (Avg from measurements)

datapoints = 100

VgsSweepFields = 'Vgs Start', 'Vgs End', 'Vgs Step'
VgsSweepParams = -10, 10, 0.2

VdsSweepFields = 'Vds Start', 'Vds End', 'Vds Step'
VdsSweepParams = 0.2, 0.2, 0.1

fields2 = 'Dielectric Thickness (nm)', 'Channel Width (um)', 'Channel Length (um)', 'Mobility', 'Phonon Energy', 'Carrier Density', 'Gate Threshold\nVoltage'
params2 = tox1, W, L, mu, Ep, N, Vgt0

default_resolution = [800, 600]
top_height = 100

class GUI:

    def __init__(self, master):
        self.root = master
        self.root.geometry(str(default_resolution[0]) + 'x' + str(default_resolution[1]))
        self.root.title('GFET Simulation')

        self.data = {}
        
        top = tk.Frame(self.root, width=default_resolution[0], height=top_height)
        left = tk.Frame(self.root, width=int(default_resolution[0]/2), height=(default_resolution[1]-top_height))
        right = tk.Frame(self.root, width=int(default_resolution[0]/2), height=(default_resolution[1]-top_height))

        
        # layout all of the main containers
        self.root.grid_rowconfigure(1, weight=1)
        self.root.grid_columnconfigure(0, weight=1)

        top.pack(fill="both")
        left.pack(side="left")
        right.pack(side="right")

        # Left Panel Tabs (settings etc)
        tab_control = ttk.Notebook(left)
        tab_control.pack()

        self.ltab1 = ttk.Frame(tab_control)
        tab_control.add(self.ltab1,text= 'Sweep Parameters')
                        
        self.ltab2 = ttk.Frame(tab_control)
        tab_control.add(self.ltab2, text= 'Device Parameters')

        self.setupParamsTab(self.ltab2, self.root)
        self.setupSweepTab(self.ltab1, self.root)

        # Right Panel Tabs (plots etc)
        tab_control2 = ttk.Notebook(right)
        tab_control2.pack()

        self.rtab1 = ttk.Frame(tab_control2)
        tab_control2.add(self.rtab1,text= 'Transfer Characteristic')
                        
        self.rtab2 = ttk.Frame(tab_control2)
        tab_control2.add(self.rtab2, text= 'I-V Characteristic')

        self.setupAxes(self.rtab1, self.rtab2)

        # Top Panel stuff (i.e. buttons)
        self.setupTopFrame(top)

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
#            lab.grid(row=index+2,column=0)
#            ent.grid(row=index+2,column=1)
            entries.append((field, ent))
        return entries

    def generateSweep(self, ent, profile, vgsModel, vdsModel):

        ent1 = self.fetch(ent)

        Vgs_start = float(ent1[0])
        Vgs_end = float(ent1[1])
        Vgs_step = float(ent1[2])
        
        Vds_start = float(ent1[3])
        Vds_end = float(ent1[4])
        Vds_step = float(ent1[5])

##        print('\nSweep Parameters:')
##        print('Vgs Range: ' + str(Vgs_start) + ' to ' + str(Vgs_end) + ', Step: ' +str(Vgs_step))
##        print('Vds Range: ' + str(Vds_start) + ' to ' + str(Vds_end) + ', Step: ' +str(Vds_step))

        # See if step correction is needed, i.e. a step must be at least 1
        if profile == "Sweep Vgs, Step Vds":
            VgsStepCorrection = 0
            VdsStepCorrection = 1
        elif profile == "Step Vgs, Sweep Vds":
            VgsStepCorrection = 1
            VdsStepCorrection = 10
            
        # Vgs Model
        if vgsModel == "Linear":
            dps = VgsStepCorrection + int(abs(Vgs_start/Vgs_step) + abs(Vgs_end/Vgs_step))
            Vgs = list(np.linspace(Vgs_start, Vgs_end, dps))
        elif vgsModel == "Dual-Linear":
            dps = VgsStepCorrection + int(abs(Vgs_start/Vgs_step) + abs(Vgs_end/Vgs_step))
            # i.e. forwards and backwards sweep
            Vgs = list(np.linspace(Vgs_start, Vgs_end, dps)) + list(np.linspace(Vgs_end, Vgs_start, dps))
        elif vgsModel == "Logarithmic":
            dps = VgsStepCorrection + int(abs(Vgs_start/Vgs_step) + abs(Vgs_end/Vgs_step))
            # i.e. forwards and backwards sweep
            Vgs = list(np.logspace(Vgs_start, Vgs_end, dps)) + list(np.linspace(Vgs_end, Vgs_start, dps))
            
        # Vds Model
        if vdsModel == "Linear":
            dps = VdsStepCorrection + int(abs(round((Vds_start - Vds_end)/Vds_step)))     
            Vds = list(np.linspace(Vds_start, Vds_end, dps))
        elif vdsModel == "Dual-Linear":
            dps = VdsStepCorrection + int(abs(round((Vds_start - Vds_end)/Vds_step)))     
            Vds = list(np.linspace(Vds_start, Vds_end, dps))
        elif vdsModel == "Logarithmic":
            dps = VdsStepCorrection + int(abs(round((Vds_start - Vds_end)/Vds_step)))     
            Vds = list(np.logspace(Vds_start, Vds_end, dps))        
        return {"Vgs": Vgs,
                "Vds": Vds}

    def loadModel(self, name, profile, vgsModel, vdsModel, ents, ents2):
        self.model = name

        # Todo: figure out loading sweeps vs getting from the GUI
        self.data.update(self.generateSweep(ents, profile, vgsModel, vdsModel))

        params = self.fetch(ents2)
        eps = self.dielecCombo

        if self.model == 'Rodriguez':
            GFET = gfet.RodriguezGFET(params, self.data["Vds"], self.data["Vgs"], eps)
            self.data.update(GFET.calculateIds())
        elif self.model == 'Thiele':
            GFET = gfet.ThieleGFET(params, self.data["Vds"], self.data["Vgs"], eps)
            self.data.update(GFET.calculateIds())
        else:
            print('ERROR: No models available.')

        if self.model != None:
            self.plotTransferChars() # maybe all curves?
            self.plotIVChars()

    def plotTransferChars(self):
        plotted = False
        
        self.ax.clear()
        self.ax.set_title('Transfer Characteristic')
        self.ax.set_ylabel(r'$I_{DS}$ (A)')
        self.ax.set_xlabel(r'$V_{GS}$ (V)')

        if "Vgs" in self.data:
            for index,entry in enumerate(self.data["Ids"]):
                self.ax.plot(self.data["Vgs"], self.data["Ids"][index])
            plotted = True

        self.ax.set_aspect(1./self.ax.get_data_ratio())
        self.ax.ticklabel_format(axis='y',style='sci', scilimits=(0,0))
        self.canvas.draw()
        
    def plotIVChars(self):
        plotted = False
        
        self.ax2.clear()
        self.ax2.set_title('I-V Characteristic')
        self.ax2.set_ylabel(r'$I_{DS}$ (A)')
        self.ax2.set_xlabel(r'$V_{DS}$ (V)')

        if "Vds" in self.data:
            for index,entry in enumerate(self.data["Vds"]):
                self.ax2.plot(self.data["Vds"], self.data["Ids"])
            plotted = True
        
        self.ax2.set_aspect(1./self.ax2.get_data_ratio())
        self.ax2.ticklabel_format(axis='y',style='sci', scilimits=(0,0))
        self.canvas2.draw() 
    
    def setupSweepTab(self, tab1, root):
        frame = tk.Frame(tab1)

        # Setup Model Selection Box
        modelLabel = tk.Label(frame, text='GFET Model')
        self.modelCombo = ttk.Combobox(frame, state='readonly')
        self.modelCombo['values'] = 'Rodriguez', 'Thiele' # Will need to auto generate, maybe can
                                                          # see if from class names in models?
        self.modelCombo.current(0) # First value in list is default
        modelLabel.grid(row=0, column=0)
        self.modelCombo.grid(row=0, column=1)
    
        frame.pack()

        # Setup Sweep Selection Box
        # Which sweep to do
        sweepOptionLabel = tk.Label(frame, text='Voltage Sweep Profile')
        self.sweepOption = ttk.Combobox(frame, state='readonly')
        self.sweepOption['values'] = 'Sweep Vgs, Step Vds', 'Step Vgs, Sweep Vds'
        self.sweepOption.current(0) # First value in list is default
        sweepOptionLabel.grid(row=1, column=0)
        self.sweepOption.grid(row=1, column=1)


        # Vgs Sweep Settings
        VgsSweepLabel = tk.Label(frame, text='Vgs Sweep Model')
        self.VgsSweepCombo = ttk.Combobox(frame, state='readonly')
        self.VgsSweepCombo['values'] = 'Linear', 'Dual-Linear', 'Logarithmic'
        self.VgsSweepCombo.current(0) # First value in list is default
        VgsSweepLabel.grid(row=2, column=0)
        self.VgsSweepCombo.grid(row=2, column=1)

        self.ents = self.makeform(tab1, VgsSweepFields, VgsSweepParams)
        root.bind('<Return>', (lambda event, e=self.ents: self.fetch(e)))

        # Vds Sweep Settings
        VdsSweepLabel = tk.Label(frame, text='Vds Sweep Model')
        self.VdsSweepCombo = ttk.Combobox(frame, state='readonly')
        self.VdsSweepCombo['values'] = 'Linear', 'Dual-Linear', 'Logarithmic'
        self.VdsSweepCombo.current(0) # First value in list is default
        VdsSweepLabel.grid(row=3, column=0)
        self.VdsSweepCombo.grid(row=3, column=1)

        self.ents = self.ents + (self.makeform(tab1, VdsSweepFields, VdsSweepParams))
        root.bind('<Return>', (lambda event, e=self.ents: self.fetch(e)))

        root.bind("<<ComboboxSelected>>")
        
    def setupTopFrame(self,top):
        io = gio.GFET_IO()
        
        b1 = tk.Button(top, text='Simulate',
                  command=(lambda e=self.ents, e2=self.ents2: self.loadModel(self.modelCombo.get(), self.sweepOption.get(), self.VgsSweepCombo.get(), self.VdsSweepCombo.get(), e, e2)))
        b1.pack(side='left')
        b2 = tk.Button(top, text='Load Sweep', command=io.loadSweep)
        b2.pack(side='left')
        b3 = tk.Menubutton(top, text='Export Data')
        b3.menu = tk.Menu(b3)
        b3["menu"] = b3.menu
        b3.menu.add_command(label="Export Transfer Characteristics", command=(lambda : io.exportTransferChars(self.data)))
        b3.menu.add_command(label="Export I-V Characteristics", command=(lambda : io.exportIVChars(self.data)))
        b3.pack(side='left')
        
    def setupParamsTab(self, tab2, root):
        frame = tk.Frame(tab2)
        dielecLabel = tk.Label(frame, text='Relative Permittivity')
        self.dielecCombo = ttk.Combobox(frame, state='readonly')

        dielecs = gio.GFET_IO.loadDielectrics()
        
        vals = []
        i = 0
        for i in range(len(dielecs)):
            vals.append(dielecs[i][0] + " (" + str(dielecs[i][1]) + ")")

        self.dielecCombo['values'] = vals
        self.dielecCombo.current(0) # First value in list is default
        
        root.bind("<<ComboboxSelected>>")
        frame.pack(side=tk.TOP, fill=tk.X, padx=5, pady=5)
        dielecLabel.pack(side=tk.LEFT)
        self.dielecCombo.pack(side=tk.RIGHT)

        subFrame = tk.Frame()
        
        self.ents2 = self.makeform(tab2, fields2, params2)
        root.bind('<Return>', (lambda event, e=self.ents2: self.fetch(e)))


    def setupAxes(self, frame1, frame2):
        fig = plt.Figure(figsize=(4,4), dpi=100)
        self.ax = fig.add_subplot(111)

        self.canvas = FigureCanvasTkAgg(fig, master=frame1)
        self.canvas.get_tk_widget().pack()

        fig2 = plt.Figure(figsize=(4,4), dpi=100)
        self.ax2 = fig2.add_subplot(111)

        self.canvas2 = FigureCanvasTkAgg(fig2, master=frame2)
        self.canvas2.get_tk_widget().pack()
        
        self.plotTransferChars()
        self.plotIVChars()
