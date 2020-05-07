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
from scipy import constants as consts

# Default Device Parameters
tox1 = 50 # TG oxide layer thickness, metres
W = 40 # Channel Width, metres
L = 120 # Channel Length, metres
mu = 7000 # Effective carrier mobility
N = 5e-17 # Dopant density, assume basically zero
Ep = 5.6e-20 # Surface phonon energy of the substrate

# Default Sweep Settings
transSweepFields = 'Vgs Start', 'Vgs End', 'Vgs Step', 'Vds Start', 'Vds End', 'Vds Step'
transSweepParams = -10, 10, 0.2,  0.2, 0.2, 0.1

ivSweepFields = 'Vgs Start', 'Vgs End', 'Vgs Step', 'Vds Start', 'Vds End', 'Vds Step'
ivSweepParams = 0, 2, 0.5, 0.1, 1, 0.01

fields2 = ('Dielectric Thickness (nm)', 'Channel Width (um)', 'Channel Length (um)', 'Mobility (cm2/V/s)',
            'Phonon Energy (J)', 'Effective Dopant\nDensity')
params2 = tox1, W, L, mu, Ep, N

default_resolution = [1024, 600]
top_height = 100

class GUI:

    def __init__(self, master):
        self.root = master
        self.root.geometry(str(default_resolution[0]) + 'x' + str(default_resolution[1]))
        self.root.title('GFET Simulator')

        self.io = gio.GFET_IO()

        self.data = {}
        
        top = tk.Frame(self.root, width=default_resolution[0], height=top_height)
        left = tk.Frame(self.root, width=int(2*default_resolution[0]/5),
                        height=(default_resolution[1]-top_height))
        right = tk.Frame(self.root, width=int(3*default_resolution[0]/5),
                         height=(default_resolution[1]-top_height))

        top.pack(fill="both")
        left.pack(side="left")
        right.pack(side="right")

        for frame in [left, right]:
            frame.pack_propagate(0)

        # Left Panel Tabs (settings etc)
        tab_control = ttk.Notebook(left, style='Custom.TNotebook')
        tab_control.pack()

        self.ltab1 = ttk.Frame(tab_control)
        tab_control.add(self.ltab1,text= 'Sweep Parameters')
                        
        self.ltab2 = ttk.Frame(tab_control)
        tab_control.add(self.ltab2, text= 'Device Parameters')

        self.setupParamsTab()
        self.setupSweepTab()

        # Right Panel Tabs (plots etc)
        tab_control2 = ttk.Notebook(right, style='Custom.TNotebook')
        tab_control2.pack()

        self.rtab1 = ttk.Frame(tab_control2)
        tab_control2.add(self.rtab1,text= 'Transfer Characteristic')
                        
        self.rtab2 = ttk.Frame(tab_control2)
        tab_control2.add(self.rtab2, text= 'I-V Characteristic')

        self.rtab3 = ttk.Frame(tab_control2)
        tab_control2.add(self.rtab3, text= 'Conductance')

        self.rtab4 = ttk.Frame(tab_control2)
        tab_control2.add(self.rtab4, text= 'Transit Frequency')

        self.setupAxes()

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
            lab = tk.Label(row, width=20, text=field, anchor='w')

            ent = tk.Entry(row, bd=1, width=25)
            ent.insert(0, str(params[index]))
            
            row.pack(side=tk.TOP, fill=tk.X, padx=5, pady=5)
            lab.pack(side=tk.LEFT)
            ent.pack(side=tk.RIGHT, expand=tk.YES, fill=tk.X)
            entries.append((field, ent))
        return entries

#**************************************************************************************************#
#                             Model & Sweep Functions                                              #
#**************************************************************************************************#

    def generateTransferSweep(self, ent, vgsModel, vdsModel):
        ent1 = self.fetch(ent)

        Vgs_start = float(ent1[0])
        Vgs_end = float(ent1[1])
        Vgs_step = float(ent1[2])
        
        Vds_start = float(ent1[3])
        Vds_end = float(ent1[4])
        Vds_step = float(ent1[5])

        VgsStepCorrection = 0
        VdsStepCorrection = 1

        retDict = self.genSweepModels(vgsModel, vdsModel, Vgs_start, Vgs_end, Vgs_step,
                                      VgsStepCorrection, Vds_start, Vds_end,
                                      Vds_step, VdsStepCorrection) 
        return retDict

    def generateIVSweep(self, ent, vgsModel, vdsModel):
        ent1 = self.fetch(ent)

        Vgs_start = float(ent1[0])
        Vgs_end = float(ent1[1])
        Vgs_step = float(ent1[2])
        Vds_start = float(ent1[3])
        Vds_end = float(ent1[4])
        Vds_step = float(ent1[5])

        VgsStepCorrection = 1
        VdsStepCorrection = 0

        retDict = self.genSweepModels(vgsModel, vdsModel, Vgs_start, Vgs_end, Vgs_step,
                                      VgsStepCorrection, Vds_start, Vds_end,
                                      Vds_step, VdsStepCorrection)            
        return retDict

    # Used by generateTransferSweep and by generateIVSweep
    def genSweepModels(self, vgsModel, vdsModel, Vgs_start, Vgs_end, Vgs_step, VgsStepCorrection,
                       Vds_start, Vds_end, Vds_step, VdsStepCorrection):

        # Vgs Model
        if vgsModel == "Linear":
            dps = VgsStepCorrection + int(abs(Vgs_start/Vgs_step) + abs(Vgs_end/Vgs_step))
            Vgs = list(np.linspace(Vgs_start, Vgs_end, dps))
        elif vgsModel == "Dual-Linear":
            dps = VgsStepCorrection + int(abs(Vgs_start/Vgs_step) + abs(Vgs_end/Vgs_step))
            # i.e. forwards and backwards sweep
            Vgs = (list(np.linspace(Vgs_start, Vgs_end, dps))
                    + list(np.linspace(Vgs_end, Vgs_start, dps)))
        elif vgsModel == "Logarithmic":
            dps = VgsStepCorrection + int(abs(Vgs_start/Vgs_step) + abs(Vgs_end/Vgs_step))
            # i.e. forwards and backwards sweep
            Vgs = (list(np.logspace(Vgs_start, Vgs_end, dps))
                    + list(np.linspace(Vgs_end, Vgs_start, dps)))
            
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
    
    def loadModel(self, name, vgsModel, vdsModel, ents, ents2, ents3):
        self.model = name

        ivSweep = self.generateIVSweep(ents2, vgsModel, vdsModel)

        if self.io.extSweep:
            transferSweep = self.io.transData
            subSweep = {"Vds": self.generateTransferSweep(ents, vgsModel, vdsModel)["Vds"]}
            transferSweep.update(subSweep)
        else: 
            transferSweep = self.generateTransferSweep(ents, vgsModel, vdsModel)
        
        self.data.update({"IVChars": ivSweep,
                          "TransChars": transferSweep})

        params = self.fetch(ents3)
        eps = self.dielecCombo

        if self.model == 'Rodriguez':
            GFET = gfet.RodriguezGFET(params, ivSweep, transferSweep, eps)
            transferChars, gm, fT  = GFET.calculateTransferChars()
            ivChars = GFET.calculateIVChars()
            self.data["IVChars"].update({"Ids": ivChars})
            self.data["TransChars"].update({"Ids": transferChars})
            self.data.update({"gm": gm,
                              "fT": fT})
        elif self.model == 'Thiele':
            GFET = gfet.ThieleGFET(params, ivSweep, transferSweep, eps)
            transferChars = GFET.calculateTransferChars()
            ivChars = GFET.calculateIVChars()
            self.data["IVChars"].update({"Ids": ivChars})
            self.data["TransChars"].update({"Ids": transferChars})
        
        if self.model != None:
            self.plotTransferChars(transferSweep["Vgs"], transferChars)
            self.plotIVChars(ivSweep["Vds"], ivChars)
            self.plotConductance(transferSweep["Vgs"], gm)
            self.plotFrequency(transferSweep["Vgs"], fT)            

    def loadSweep(self):
        self.io.loadSweep()

        for entry in self.ents:
            entry[1].config(state="disabled")

#**************************************************************************************************#
#                             Plotting Functions                                                   #
#**************************************************************************************************#

    # Plot the transfer characteristics
    def plotTransferChars(self, Vgs, Ids):
        plotted = False
        
        self.ax.clear()
        self.ax.set_title('Drain Current  vs Vgs')
        self.ax.set_ylabel(r'$I_{DS}$ (A)')
        self.ax.set_xlabel(r'$V_{GS}$ (V)')

        if Vgs:
            for index,entry in enumerate(Ids):
                 self.ax.plot(Vgs, Ids[index])
            plotted = True

        self.ax.set_aspect(1./self.ax.get_data_ratio())
        self.ax.ticklabel_format(axis='y',style='sci', scilimits=(0,0))
        self.canvas.draw()

    # Plot the IV Characteristics
    def plotIVChars(self, Vds, Ids):
        plotted = False
        
        self.ax2.clear()
        self.ax2.set_title('Drain Current vs Vds')
        self.ax2.set_ylabel(r'$I_{DS}$ (A)')
        self.ax2.set_xlabel(r'$V_{DS}$ (V)')

        if Vds:
            for index,entry in enumerate(Ids):
                self.ax2.plot(Vds, Ids[index])
            plotted = True
        
        self.ax2.set_aspect(1./self.ax2.get_data_ratio())
        self.ax2.ticklabel_format(axis='y',style='sci', scilimits=(0,0))
        self.canvas2.draw()

    # Plot the Conductance Characteristics
    def plotConductance(self, Vgs, gm):
        plotted = False
        
        self.ax3.clear()
        self.ax3.set_title('Conductance vs Vgs')
        self.ax3.set_ylabel(r'$g_m$ (S)')
        self.ax3.set_xlabel(r'$V_{GS}$ (V)')

        if Vgs:
            for entry in gm:
                self.ax3.plot(Vgs, entry)
                plotted = True
            
        self.ax3.set_aspect(1./self.ax3.get_data_ratio())
        self.ax3.ticklabel_format(axis='y',style='sci', scilimits=(0,0))
        self.canvas3.draw()

    def plotFrequency(self, Vgs, fT):
        plotted = False

        self.ax4.clear()
        self.ax4.set_title('Transit Frequency vs Vgs')
        self.ax4.set_ylabel(r'$f_T$ (Hz)')
        self.ax4.set_xlabel(r'$V_{GS}$ (V)')

        if Vgs:
            for entry in fT:
                self.ax4.plot(Vgs, entry)
                plotted = True
            
        self.ax4.set_aspect(1./self.ax4.get_data_ratio())
        self.ax4.ticklabel_format(axis='y',style='sci', scilimits=(0,0))
        self.canvas4.draw()
    
    def restoreDefaultSettings(self, notebook):
        current = notebook.index("current")
                
        if current == 0: # Trans sweep
            for index,entry in enumerate(self.ents):
                entry[1].config(state="normal")
                entry[1].delete(0, len(entry[1].get()))
                entry[1].insert(0, transSweepParams[index])
        elif current == 1: # IV sweep
            for index,entry in enumerate(self.ents2):
                entry[1].config(state="normal")
                entry[1].select_range(0, len(self.fetch(entry[1])))
                entry[1].select_clear()
                entry[1].insert(0, ivSweepParams[index])

#**************************************************************************************************#
#                             GUI Setup Functions                                                  #
#**************************************************************************************************#
                
    # Tab for voltage sweep settings
    def setupSweepTab(self):     
        top = tk.Frame(self.ltab1, height=(self.ltab1.winfo_height()/2))
        bottom = tk.Frame(self.ltab1, height=(self.ltab1.winfo_height()/2))

        # Setup Model Selection Box
        modelLabel = tk.Label(top, text='GFET Model')
        self.modelCombo = ttk.Combobox(top, state='readonly')
        self.modelCombo['values'] = 'Rodriguez', 'Thiele' 
        self.modelCombo.current(0) # First value in list is default
        modelLabel.grid(row=0, column=0)
        self.modelCombo.grid(row=0, column=1)
    
        top.pack()
        bottom.pack()

        # Vgs Sweep Settings
        VgsSweepLabel = tk.Label(top, text='Vgs Sweep Model')
        self.VgsSweepCombo = ttk.Combobox(top, state='readonly')
        self.VgsSweepCombo['values'] = 'Linear', 'Dual-Linear', 'Logarithmic'
        self.VgsSweepCombo.current(0) # First value in list is default
        VgsSweepLabel.grid(row=2, column=0)
        self.VgsSweepCombo.grid(row=2, column=1)

        # Vds Sweep Settings
        VdsSweepLabel = tk.Label(top, text='Vds Sweep Model')
        self.VdsSweepCombo = ttk.Combobox(top, state='readonly')
        self.VdsSweepCombo['values'] = 'Linear', 'Dual-Linear', 'Logarithmic'
        self.VdsSweepCombo.current(0) # First value in list is default
        VdsSweepLabel.grid(row=3, column=0)
        self.VdsSweepCombo.grid(row=3, column=1)

        self.root.bind("<<ComboboxSelected>>")

        # Tabs for transfer chars sweep and for IV chars sweep
        tab_control = ttk.Notebook(bottom)
        tab_control.pack()

        stab1 = ttk.Frame(tab_control)
        tab_control.add(stab1, text= 'Transfer Characteristics')
                        
        tab2 = ttk.Frame(tab_control)
        tab_control.add(tab2, text= 'I-V Characteristics')

        self.ents = self.makeform(stab1, transSweepFields, transSweepParams)
        self.ents2 = self.makeform(tab2, ivSweepFields, ivSweepParams)

        b1 = tk.Button(bottom, text='Reset Sweep',
                       command=(lambda : self.restoreDefaultSettings(tab_control)))
        b1.pack()
        
        self.root.bind('<Return>', (lambda event, e=self.ents: self.fetch(e)))

    def setupTopFrame(self,top):
        b1 = tk.Button(top, text='Simulate',
                  command=(lambda e=self.ents, e2=self.ents2, e3 = self.ents3:
                           self.loadModel(self.modelCombo.get(), self.VgsSweepCombo.get(),
                           self.VdsSweepCombo.get(), e, e2, e3)))
        b1.pack(side='left')
        b2 = tk.Button(top, text='Load Sweep', command=self.loadSweep)
        b2.pack(side='left')
        b3 = tk.Menubutton(top, text='Export Data')
        b3.menu = tk.Menu(b3)
        b3["menu"] = b3.menu
        b3.menu.add_command(label="Export Transfer Characteristics",
                            command=(lambda : self.io.exportTransferChars(self.data)))
        b3.menu.add_command(label="Export I-V Characteristics",
                            command=(lambda : self.io.exportIVChars(self.data)))
        b3.pack(side='left')
        
    def setupParamsTab(self):        
        frame = tk.Frame(self.ltab2)
        dielecLabel = tk.Label(frame, text='Relative Permittivity')
        self.dielecCombo = ttk.Combobox(frame, state='readonly')

        dielecs = self.io.loadDielectrics()
        
        vals = []
        i = 0
        for i in range(len(dielecs)):
            vals.append(dielecs[i][0] + " (" + str(dielecs[i][1]) + ")")

        self.dielecCombo['values'] = vals
        self.dielecCombo.current(0) # First value in list is default
        
        self.root.bind("<<ComboboxSelected>>")
        frame.pack(side=tk.TOP, fill=tk.X, padx=5, pady=5)
        dielecLabel.pack(side=tk.LEFT)
        self.dielecCombo.pack(side=tk.RIGHT)

        subFrame = tk.Frame()
        
        self.ents3 = self.makeform(self.ltab2, fields2, params2)
        self.root.bind('<Return>', (lambda event, e=self.ents3: self.fetch(e)))


    def setupAxes(self):
        
        # Transfer Characteristic Plot
        fig = plt.Figure(figsize=(5,5), dpi=100)
        self.ax = fig.add_subplot(111)

        self.canvas = FigureCanvasTkAgg(fig, master=self.rtab1)
        self.canvas.get_tk_widget().pack()

        # I-V Characteristic Plot
        fig2 = plt.Figure(figsize=(5,5), dpi=100)
        self.ax2 = fig2.add_subplot(111)

        self.canvas2 = FigureCanvasTkAgg(fig2, master=self.rtab2)
        self.canvas2.get_tk_widget().pack()

        # Conductance Plot
        fig3 = plt.Figure(figsize=(5,5), dpi=100)
        self.ax3 = fig3.add_subplot(111)

        self.canvas3 = FigureCanvasTkAgg(fig3, master=self.rtab3)
        self.canvas3.get_tk_widget().pack()

        # Transit Frequency Plot
        fig4 = plt.Figure(figsize=(5,5), dpi=100)
        self.ax4 = fig4.add_subplot(111)

        self.canvas4 = FigureCanvasTkAgg(fig4, master=self.rtab4)
        self.canvas4.get_tk_widget().pack()

        # Initialise empty plot
        self.plotTransferChars(None, None)
        self.plotIVChars(None, None)
        self.plotConductance(None, None)
        self.plotFrequency(None, None)
