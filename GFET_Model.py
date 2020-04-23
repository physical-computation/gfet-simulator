# Program to model electrical properties of Graphene Field-Effect Transistors (GFETs).
# Based on S. Rodriguez et al., "A Comprehensive Graphene FET Model for Circuit Design,"
# in IEEE Transactions on Electron Devices, vol. 61, no. 4, pp. 1199-1206, April 2014. (doi.org:/10.1109/TED.2014.2302372)

import numpy as np

from numpy.lib.scimath import sqrt as csqrt

from scipy import constants as consts # For Physical Constants

import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import (
    FigureCanvasTkAgg, NavigationToolbar2Tk)
# Implement the default Matplotlib key bindings.
from matplotlib.backend_bases import key_press_handler
from matplotlib.figure import Figure
import matplotlib.animation as anim

import tkinter as tk
from tkinter import ttk

#******************************************************************************#
#                                   Main Program                               #
#******************************************************************************#


# Device Parameters
er = 9.34; # Oxide permittivity, Alumina
tox1 = 50 # TG oxide layer thickness, metres
W = 40 # Channel Width, metres
L = 120 # Channel Length, metres
mu = 7000 # Effective carrier mobility
N = 0.5*10**(16) # carrier density
Ep = 0.056*10**(-19) # Surface phonon energy of the substrate
Vgt0 = -0.8229 # TG Voltage at Dirac pt, volts. (Avg from measurements)

# Initial Calculations
Ct = 3.6*10**(-3)
w = Ep/consts.hbar

fields1 = '-Vgs', '+Vgs', 'Vds', 'Datapoints'
params1 = -10, 10, 0.2, 100

fields2 = 'Dielectric Permittivity', 'Dielectric Thickness (nm)', 'Channel Width (um)', 'Channel Length (um)', 'Mobility'
params2 = er, tox1, W, L, mu

V = []
I = []

def fetch(entries):
    data = []
    for entry in entries:
        text  = entry[1].get()
        data.append(text)

    return data

def makeform(root, fields, params):
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

def plotCurve(V, I, ax, canvas):
    ax.clear()
    ax.set_title('GFET Transfer Characteristic')
    ax.plot(V, I)
    canvas.draw()

def init_func(i):
    ax.clear()
    ax.set_title('GFET Transfer Characteristic')

def calculateIds(entries, ax, canvas):
    
    data = fetch(entries)
    start_voltage = float(data[0])
    end_voltage = float(data[1])
    Vds = float(data[2])
    datapoints = int(data[3])
    
    Vgt = np.linspace(start_voltage, end_voltage, datapoints)
    Id = []

    # Convert values to correct units
    Wr = W*10**(-6)
    Lr = L*10**(-6)
    toxr = tox1*10**(-9)
    
    for i in range(datapoints):
        Veff = Vgt[i] + Vgt0
        Id.append(abs((mu*Wr*Ct*(Veff-0.5*Vds))/(Lr/Vds + (mu/w)*(csqrt(consts.pi*Ct/consts.elementary_charge))*(csqrt(Veff-0.5*Vds)))))

    plotCurve(Vgt, Id, ax, canvas)
    
if __name__ == '__main__':

    root = tk.Tk()
    root.geometry('800x600')
    root.title('GFET Simulation')

    top = tk.Frame(root, width=640, height=50)
    left = tk.Frame(root, width=220, height=430)
    right = tk.Frame(root, width=420, height=430)

    # layout all of the main containers
    root.grid_rowconfigure(1, weight=1)
    root.grid_columnconfigure(0, weight=1)

    top.pack(side="top")
    left.pack(side="left")
    right.pack(side="right")

    # Setup Left Frame Tabs    

    tab_control = ttk.Notebook(left)
    tab_control.pack() # expand=1, fill='both')

    tab1 = ttk.Frame(tab_control)
    tab_control.add(tab1,text= 'Voltage Sweep')
                    
    tab2 = ttk.Frame(tab_control)
    tab_control.add(tab2, text= 'Device Parameters')

    # Additional Widgets    
    ents = makeform(tab1, fields1, params1)
    root.bind('<Return>', (lambda event, e=ents: fetch(e)))

    ents2 = makeform(tab2, fields2, params2)
    root.bind('<Return>', (lambda event, e=ents2: fetch(e)))
    
    fig = plt.Figure(figsize=(5,4), dpi=100)
    ax = fig.add_subplot(111)
    ax.set_title('GFET Transfer Characteristic')
      
    b1 = tk.Button(tab1, text='Simulate',
                  command=(lambda e=ents: calculateIds(e, ax, canvas)))
    b1.pack(side=tk.LEFT, padx=5, pady=5)
    b2 = tk.Button(tab1, text='Quit', command=root.quit)
    b2.pack(side=tk.LEFT, padx=5, pady=5)

    canvas = FigureCanvasTkAgg(fig, master=right)
    canvas.get_tk_widget().pack()
    canvas.draw()

    root.mainloop()
    
