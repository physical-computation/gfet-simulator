# Class to setup and manage GUI components for GFET simulator program.

import GFET_IO as gio
import models as gfet

import tkinter as tk
from tkinter import ttk, filedialog

import re

import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import (
    FigureCanvasTkAgg, NavigationToolbar2Tk)
from matplotlib.backend_bases import key_press_handler
from matplotlib.figure import Figure
from matplotlib.widgets import Cursor

import numpy as np
from scipy import constants as consts

# Default Device Parameters
tox1 = 50 # TG oxide layer thickness, nm
tox2 = 50 # BG oxide layer thickness, nm
W = 40 # Channel Width, um
L = 120 # Channel Length, um
mu = 7000 # Effective carrier mobility
N = 5e-17 # Dopant density, assume basically zero
Ep = 5.6e-20 # Surface phonon energy of the substrate
T = 298 # Operating temperature, default is room temp (in Kelvin)

# Default Sweep Settings
transSweepFields = 'Vtg Start', 'Vtg End', 'Vtg Step', 'Vbg Start', 'Vbg End', 'Vbg Step', 'Vds Start', 'Vds End', 'Vds Step'
transSweepParams = -10, 10, 0.2, 0, 0, 1,  0.2, 0.2, 0.1

ivSweepFields = 'Vtg Start', 'Vtg End', 'Vtg Step', 'Vbg Start', 'Vbg End', 'Vbg Step', 'Vds Start', 'Vds End', 'Vds Step'
ivSweepParams = 0, 2, 0.5, 0, 0, 1, 0, 1, 0.01

fields2 = ('Top Dielectric\nThickness (nm)', 'Bottom Dielectric\nThickness (nm)', 'Channel Width (um)', 'Channel Length (um)', 'Mobility (cm2/V/s)',
            'Phonon Energy (J)', 'Effective Dopant\nDensity', 'Operating Temperature (K)')
params2 = tox1, tox2, W, L, mu, Ep, N, T

default_resolution = [1024, 600]
top_height = 100
scatter_size = 10

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
        tab_control2.add(self.rtab3, text= 'Transconductance')

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

    # See if the entry can be converted to a float or not. If not, it's not valid input
    def validate_entry(self, content, newcont):
        try:
            float(content)
            return True
        except ValueError:
            return False

    def makeform(self, root, fields, params):
        entries = []

        validate = (root.register(self.validate_entry), "%P", '%s')

        for index,field in enumerate(fields):
            row = tk.Frame(root)
            lab = tk.Label(row, width=20, text=field, anchor='w')
            
            ent = tk.Entry(row, bd=1, width=25, validate="key", validatecommand=validate)
            ent.insert(0, str(params[index]))

            row.pack(side=tk.TOP, fill=tk.X, padx=5, pady=5)
            lab.pack(side=tk.LEFT)
            ent.pack(side=tk.RIGHT, expand=tk.YES, fill=tk.X)
            entries.append((field, ent))
        return entries

#**************************************************************************************************#
#                             Model & Sweep Functions                                              #
#**************************************************************************************************#

    def generateTransferSweep(self, ent, vtgModel, vbgModel, vdsModel):
        ent1 = self.fetch(ent)

        Vtg_start = float(ent1[0])
        Vtg_end = float(ent1[1])
        Vtg_step = float(ent1[2])
        
        Vbg_start = float(ent1[3])
        Vbg_end = float(ent1[4])
        Vbg_step = float(ent1[5])
        
        Vds_start = float(ent1[6])
        Vds_end = float(ent1[7])
        Vds_step = float(ent1[8])

        VtgStepCorrection = 0
        VbgStepCorrection = 1
        VdsStepCorrection = 1

        retDict = self.genSweepModels(vtgModel, vbgModel, vdsModel, Vtg_start, Vtg_end, Vtg_step,
                                      Vbg_start, Vbg_end, Vbg_step,
                                      VtgStepCorrection, VbgStepCorrection, Vds_start, Vds_end,
                                      Vds_step, VdsStepCorrection) 
        return retDict

    def generateIVSweep(self, ent, vtgModel, vbgModel, vdsModel):
        ent1 = self.fetch(ent)

        Vtg_start = float(ent1[0])
        Vtg_end = float(ent1[1])
        Vtg_step = float(ent1[2])
        Vbg_start = float(ent1[3])
        Vbg_end = float(ent1[4])
        Vbg_step = float(ent1[5])
        
        Vds_start = float(ent1[6])
        Vds_end = float(ent1[7])
        Vds_step = float(ent1[8])

        VtgStepCorrection = 1
        VbgStepCorrection = 0
        VdsStepCorrection = 0

        retDict = self.genSweepModels(vtgModel, vbgModel, vdsModel, Vtg_start, Vtg_end, Vtg_step,
                                      Vbg_start, Vbg_end, Vbg_step,
                                      VtgStepCorrection, VbgStepCorrection, Vds_start, Vds_end,
                                      Vds_step, VdsStepCorrection)            
        return retDict

    # Used by generateTransferSweep and by generateIVSweep
    def genSweepModels(self, vtgModel, vbgModel, vdsModel, Vtg_start, Vtg_end, Vtg_step,
                       Vbg_start, Vbg_end, Vbg_step, VtgStepCorrection, VbgStepCorrection,
                       Vds_start, Vds_end, Vds_step, VdsStepCorrection):

        # Vtg Model
        if vtgModel == "Linear":
            dps = VtgStepCorrection + int(abs(Vtg_start/Vtg_step) + abs(Vtg_end/Vtg_step))
            Vtg = list(np.linspace(Vtg_start, Vtg_end, dps))
        elif vtgModel == "Dual-Linear":
            dps = VtgStepCorrection + int(abs(Vtg_start/Vtg_step) + abs(Vtg_end/Vtg_step))
            # i.e. forwards and backwards sweep
            Vtg = (list(np.linspace(Vtg_start, Vtg_end, dps))
                    + list(np.linspace(Vtg_end, Vtg_step, dps)))
        elif vtgModel == "Logarithmic":
            dps = VtgStepCorrection + int(abs(Vtg_start/Vtg_step) + abs(Vtg_end/Vtg_step))
            # i.e. forwards and backwards sweep
            Vtg = (list(np.logspace(Vtg_start, Vtg_end, dps))
                    + list(np.linspace(Vtg_end, Vtg_start, dps)))

        # Vbg Model
        if vbgModel == "Linear":
            dps = VbgStepCorrection + int(abs(Vbg_start/Vbg_step) + abs(Vbg_end/Vbg_step))
            Vbg = list(np.linspace(Vbg_start, Vbg_end, dps))
        elif vbgModel == "Dual-Linear":
            dps = VbgStepCorrection + int(abs(Vbg_start/Vbg_step) + abs(Vbg_end/Vbg_step))
            # i.e. forwards and backwards sweep
            Vbg = (list(np.linspace(Vbg_start, Vbg_end, dps))
                    + list(np.linspace(Vbg_end, Vbg_step, dps)))
        elif vbgModel == "Logarithmic":
            dps = VbgStepCorrection + int(abs(Vbg_start/Vbg_step) + abs(Vbg_end/Vbg_step))
            # i.e. forwards and backwards sweep
            Vbg = (list(np.logspace(Vbg_start, Vbg_end, dps))
                    + list(np.linspace(Vbg_end, Vbg_start, dps)))
        
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


        return {"Vtg": Vtg,
                "Vbg": Vbg,
                "Vds": Vds}
    
    def loadModel(self, name, vtgModel, vbgModel, vdsModel, ents, ents2, ents3):
        self.model = name

        ivSweep = self.generateIVSweep(ents2, vtgModel, vbgModel, vdsModel)

        if self.io.extSweep:
            transferSweep = self.io.transData
            subSweep = {"Vds": self.generateTransferSweep(ents, vtgModel, vbgModel, vdsModel)["Vds"]}
            transferSweep.update(subSweep)
        else: 
            transferSweep = self.generateTransferSweep(ents, vtgModel, vbgModel, vdsModel)
        
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
            ivChars = GFET.calculateIVChars()
            transferChars, gm, fT = GFET.calculateTransferChars()
            self.data["IVChars"].update({"Ids": ivChars})
            self.data["TransChars"].update({"Ids": transferChars})
        elif self.model == 'Hu':
            GFET = gfet.HuGFET(params, ivSweep, transferSweep, eps)
            transferChars, gm, fT = GFET.calculateTransferChars()
            ivChars = GFET.calculateIVChars()
            self.data["IVChars"].update({"Ids": ivChars})
            self.data["TransChars"].update({"Ids": transferChars})

        # Plot data
        if self.model != None:
            self.plotTransferChars(transferSweep["Vtg"], transferChars)
            self.plotIVChars(ivSweep["Vds"], ivChars)
            self.plotConductance(transferSweep["Vtg"], gm)
            self.plotFrequency(transferSweep["Vtg"], fT)            

    def loadSweep(self):
        self.io.loadSweep()

        for entry in self.ents:
            entry[1].config(state="disabled")

#**************************************************************************************************#
#                             Plotting Functions                                                   #
#**************************************************************************************************#

    # Plot the transfer characteristics
    def plotTransferChars(self, Vtg, Ids):
        plotted = False
        
        self.ax.clear()
        self.ax.set_title('Drain Current  vs Vtg')
        self.ax.set_ylabel(r'$I_{DS}$ (A)')
        self.ax.set_xlabel(r'$V_{TG}$ (V)')
        self.ax.set_ylim(bottom=0, top=0.01)
        
        if Vtg:
            self.ax.autoscale(enable=True)
            maxId = 0
            minId = 0
            for index,entry in enumerate(Ids):
                self.ax.scatter(Vtg, Ids[index], s=scatter_size)
                localmax = max(entry)
                localmin = min(entry)
                if localmax >= maxId:
                    maxId = localmax
                if localmin <= minId:
                    minId = localmin
            plotted = True
        else:
            self.ax.fill()
       
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
        self.ax2.set_ylim(bottom=0, top=0.01)
        
        if Vds:
            for index, entry in enumerate(Ids):
                for entry2 in Ids[index]:
                    if entry2 < 0:
                        print(entry2)
            self.ax2.autoscale(enable=True)
            maxId = 0
            minId = 0
            for index,entry in enumerate(Ids):
                self.ax2.scatter(Vds, Ids[index], s=scatter_size)
                localmax = max(entry)
                localmin = min(entry)
                if localmax >= maxId:
                    maxId = localmax
                if localmin <= minId:
                    minId = localmin
            plotted = True
        
        self.ax2.set_aspect(1./self.ax2.get_data_ratio())
        self.ax2.ticklabel_format(axis='y',style='sci', scilimits=(0,0))
        self.canvas2.draw()

    # Plot the Conductance Characteristics
    def plotConductance(self, Vtg, gm):
        plotted = False
        
        self.ax3.clear()
        self.ax3.set_title('Transconductance vs Vtg')
        self.ax3.set_ylabel(r'$g_m$ (S)')
        self.ax3.set_xlabel(r'$V_{TG}$ (V)')

        if Vtg:
            for entry in gm:
                self.ax3.scatter(Vtg, entry, s=scatter_size)
                plotted = True
            
        self.ax3.set_aspect(1./self.ax3.get_data_ratio())
        self.ax3.ticklabel_format(axis='y',style='sci', scilimits=(0,0))
        self.canvas3.draw()

    def plotFrequency(self, Vtg, fT):
        plotted = False

        self.ax4.clear()
        self.ax4.set_title('Transit Frequency vs Vtg')
        self.ax4.set_ylabel(r'$f_T$ (Hz)')
        self.ax4.set_xlabel(r'$V_{TG}$ (V)')

        if Vtg:
            for entry in fT:
                self.ax4.scatter(Vtg, entry, s=scatter_size)
                plotted = True
            
        self.ax4.set_aspect(1./self.ax4.get_data_ratio())
        self.ax4.ticklabel_format(axis='y',style='sci', scilimits=(0,0))
        self.canvas4.draw()

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

#**************************************************************************************************#
#                             GUI Setup Functions                                                  #
#**************************************************************************************************#
    
    # Tab for voltage sweep settings
    def setupSweepTab(self):
        # Disable backgate fields if no backgate
        def modelBox(*args):
            if self.modelCombo.get() == 'Rodriguez':
                for entry in self.ents:
                    if "bg" in entry[0]:
                        entry[1].config(state="disabled")
            else:
                for entry in self.ents:
                    entry[1].config(state="normal")                
        
        top = tk.Frame(self.ltab1, height=(self.ltab1.winfo_height()/2))
        bottom = tk.Frame(self.ltab1, height=(self.ltab1.winfo_height()/2))

        # Setup Model Selection Box
        modelLabel = tk.Label(top, text='GFET Model')
        self.modelCombo = ttk.Combobox(top, state='readonly')
        self.modelCombo['values'] = 'Rodriguez', 'Thiele', 'Hu'
        self.modelCombo.current(0) # First value in list is default
        modelLabel.grid(row=0, column=0)
        self.modelCombo.grid(row=0, column=1)
    
        top.pack()
        
        # Vtg Sweep Settings
        VtgSweepLabel = tk.Label(top, text='Vtg Sweep Model')
        self.VtgSweepCombo = ttk.Combobox(top, state='readonly')
        self.VtgSweepCombo['values'] = 'Linear', 'Dual-Linear', 'Logarithmic'
        self.VtgSweepCombo.current(0) # First value in list is default
        VtgSweepLabel.grid(row=2, column=0)
        self.VtgSweepCombo.grid(row=2, column=1)

        # Vbg Sweep Settings
        VbgSweepLabel = tk.Label(top, text='Vbg Sweep Model')
        self.VbgSweepCombo = ttk.Combobox(top, state='readonly')
        self.VbgSweepCombo['values'] = 'Linear', 'Dual-Linear', 'Logarithmic'
        self.VbgSweepCombo.current(0) # First value in list is default
        VbgSweepLabel.grid(row=3, column=0)
        self.VbgSweepCombo.grid(row=3, column=1)

        # Vds Sweep Settings
        VdsSweepLabel = tk.Label(top, text='Vds Sweep Model')
        self.VdsSweepCombo = ttk.Combobox(top, state='readonly')
        self.VdsSweepCombo['values'] = 'Linear', 'Dual-Linear', 'Logarithmic'
        self.VdsSweepCombo.current(0) # First value in list is default
        VdsSweepLabel.grid(row=4, column=0)
        self.VdsSweepCombo.grid(row=4, column=1)

        self.root.bind("<<ComboboxSelected>>", modelBox)

        # Tabs for transfer chars sweep and for IV chars sweep
        tab_control = ttk.Notebook(bottom)
        tab_control.pack()

        stab1 = ttk.Frame(tab_control)
        tab_control.add(stab1, text= 'Transfer Characteristics')
                        
        tab2 = ttk.Frame(tab_control)
        tab_control.add(tab2, text= 'I-V Characteristics')

        # Make transfer tab scrollable
        canvas = tk.Canvas(stab1)
        scrollbar = ttk.Scrollbar(stab1, orient='vertical', command=canvas.yview)
        scrollable_frame = ttk.Frame(canvas)

        scrollable_frame.bind("<Configure>", lambda e: canvas.configure(scrollregion=canvas.bbox('all')))
        canvas.create_window((0,0), window=scrollable_frame, anchor='nw')
        canvas.configure(yscrollcommand=scrollbar.set)

        canvas2 = tk.Canvas(tab2)
        scrollbar2 = ttk.Scrollbar(tab2, orient='vertical', command=canvas2.yview)
        scrollable_frame2 = ttk.Frame(canvas2)

        scrollable_frame2.bind("<Configure>", lambda e: canvas2.configure(scrollregion=canvas2.bbox('all')))
        canvas2.create_window((0,0), window=scrollable_frame2, anchor='nw')
        canvas2.configure(yscrollcommand=scrollbar2.set)

        bottom.pack()
        canvas.pack(side='left', fill='both', expand=True)
        scrollbar.pack(side='right', fill='y')
        canvas2.pack(side='left', fill='both', expand=True)
        scrollbar2.pack(side='right', fill='y')

        self.ents = self.makeform(scrollable_frame, transSweepFields, transSweepParams)
        self.ents2 = self.makeform(scrollable_frame2, ivSweepFields, ivSweepParams)

        b1 = tk.Button(bottom, text='Reset Sweep',
                       command=(lambda : self.restoreDefaultSettings(tab_control)))
        b1.pack()
        
        self.root.bind('<Return>', (lambda event, e=self.ents: self.fetch(e)))

    # Setup the device parameters tab
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

    # Setup the top frame
    def setupTopFrame(self,top):
        b1 = tk.Button(top, text='Simulate',
                  command=(lambda e=self.ents, e2=self.ents2, e3 = self.ents3:
                           self.loadModel(self.modelCombo.get(), self.VtgSweepCombo.get(),
                           self.VtgSweepCombo.get(), self.VdsSweepCombo.get(), e, e2, e3)))
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


    # If external sweep loaded, allows normal sweep to be run
    def restoreDefaultSettings(self, notebook):
        current = notebook.index("current")
        self.io.extSweep = False
                
        if current == 0: # Trans sweep
            for index,entry in enumerate(self.ents):
                entry[1].config(state="normal", validate="key")
                entry[1].delete(0, len(entry[1].get()))
                entry[1].insert(0, transSweepParams[index])
        elif current == 1: # IV sweep
            for index,entry in enumerate(self.ents2):
                entry[1].config(state="normal", validate="key")
                entry[1].select_range(0, len(self.fetch(entry[1])))
                entry[1].select_clear()
                entry[1].insert(0, ivSweepParams[index])
