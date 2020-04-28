# Program to model electrical properties of Graphene Field-Effect Transistors (GFETs).

import numpy as np

from numpy.lib.scimath import sqrt as csqrt

from scipy import constants as consts # For Physical Constants

import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import (
    FigureCanvasTkAgg, NavigationToolbar2Tk)
# Implement the default Matplotlib key bindings.
from matplotlib.backend_bases import key_press_handler
from matplotlib.figure import Figure
from matplotlib.widgets import Cursor

from mpldatacursor import datacursor

import tkinter as tk
from tkinter import ttk, filedialog

#******************************************************************************#
#                                   Main Program                               #
#******************************************************************************#

# Software Parameters
default_resolution = [800, 600]

# Device Parameters
tox1 = 50 # TG oxide layer thickness, metres
W = 40 # Channel Width, metres
L = 120 # Channel Length, metres
mu = 7000 # Effective carrier mobility
N = 0.5*10**(16) # carrier density
Ep = 0.56*10**(-19) # Surface phonon energy of the substrate
Vgt0 = -0.8229 # TG Voltage at Dirac pt, volts. (Avg from measurements)

V = []
I = []
extSweep = False
datapoints = 100

fields1 = '-Vgs', '+Vgs', 'Vds', 'Datapoints'
params1 = -10, 10, 0.2, datapoints

fields2 = 'Dielectric Thickness (nm)', 'Channel Width (um)', 'Channel Length (um)', 'Mobility', 'Phonon Energy', 'Carrier Density', 'Gate Threshold\nVoltage'
params2 = tox1, W, L, mu, Ep, N, Vgt0

# Load Values (dielectric materials etc. from txt file
def loadDielectrics():
    dielectrics = np.loadtxt('Dielectrics.txt', dtype=np.dtype('O'), delimiter=',', skiprows=1)
    return dielectrics

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

def exportData(V, I):
    f = filedialog.asksaveasfile(mode='w', defaultextension=".csv")
    if f is None:
        return

    # Group voltages and currents into pairs
    dataToSave = ''.join([str(i)+","+str(j)+"\n" for i,j in zip(V,I)])
    
    f.write(dataToSave)
    f.close()

def loadSweep():
    f = filedialog.askopenfile(mode='r', filetypes=[('CSV Files', '*.csv')])

    if f is not None:
        content = f.read().split('\n')
        
        global V
        global extSweep
        global datapoints
        V = [float(item) for item in content]
        print(V)
        extSweep = True
        datapoints = len(V)

def plotCurve(V, I, ax, canvas):
    ax.clear()
    ax.set_title('GFET Transfer Characteristic')
    ax.set_ylabel(r'$I_{DS}$ (A)')
    ax.set_xlabel(r'$V_{GS}$ (V)')
    lines = ax.plot(V, I)
    ax.set_aspect(1./ax.get_data_ratio())
    ax.ticklabel_format(axis='both',style='sci', scilimits=(0,0))    
#    datacursor(lines, display='single')
    canvas.draw()

def setV(v):
    global V
    V=v

def setI(i):
    global I
    I=i

def calculateIds(entries, entries2, mod, eps, ax, canvas):

    # Get data from the entry box and convert to right form for calculations   
    data = fetch(entries)
    data2 = fetch(entries2)
    
    start_voltage = float(data[0])
    end_voltage = float(data[1])
    Vds = float(data[2])

    tox = float(data2[0])*10**(-9)
    W = float(data2[1])*10**(-6)
    L = float(data2[2])*10**(-6)
    mu = float(data2[3])
    Ep = float(data2[4])
    N = float(data2[5])
    
    er = float(eps.get().split("(")[1].replace(")",""))

    Ct = er*consts.epsilon_0/tox
#    w = Ep/consts.hbar
    w = (2.24*10**(13))/consts.pi
    
    if extSweep == False:
        dps = int(data[3])
        Vgt = np.linspace(start_voltage, end_voltage, dps)
    else:
        global V
        global datapoints
        Vgt = V
        dps = datapoints
            
    Id = []

    model = mod.get()
    
    if model == 'Rodriguez':
        for i in range(dps):
            Veff = Vgt[i] + Vgt0
            Id.append(abs((mu*W*Ct*(Veff-0.5*Vds))/(L/Vds + (mu/w)*(csqrt(consts.pi*Ct/consts.elementary_charge))*(csqrt(Veff-0.5*Vds)))))
    elif model == 'Thiele':
        for i in range(dps):
            Vch = Vgt[i]-Vgt0
            num = (mu*W*(Vch**2)*Vds*consts.elementary_charge**3)/(consts.pi*(consts.hbar*10**6)**2)
            den = L - (mu*Vds/w)*(consts.pi*Vch*(2*Vch*consts.elementary_charge**2)/(consts.pi*(consts.hbar*10**6)**2))**0.5
            Id.append(abs(num/den))
    
    setV(Vgt)
    setI(Id)

    plotCurve(Vgt, Id, ax, canvas)
    
if __name__ == '__main__':

    dielecs = loadDielectrics()

    root = tk.Tk()
    root.geometry(str(default_resolution[0]) + 'x' + str(default_resolution[1]))
    root.title('GFET Simulation')

#    top_frame_height = int((default_resolution[1]/6))

#    top = tk.Frame(root, width=default_resolution[0], height=top_frame_height)
    left = tk.Frame(root, width=int(default_resolution[0]/2), height=default_resolution[1])
    right = tk.Frame(root, width=int(default_resolution[0]/2), height=default_resolution[1])

    # layout all of the main containers
    root.grid_rowconfigure(1, weight=1)
    root.grid_columnconfigure(0, weight=1)

#    top.pack(side="top")
    left.pack(side="left")
    right.pack(side="right")

    # Setup Left Frame Tabs    

    tab_control = ttk.Notebook(left)
    tab_control.pack() # expand=1, fill='both')

    tab1 = ttk.Frame(tab_control)
    tab_control.add(tab1,text= 'Sweep Parameters')
                    
    tab2 = ttk.Frame(tab_control)
    tab_control.add(tab2, text= 'Device Parameters')

    # Additional Widgets
    ents = makeform(tab1, fields1, params1)
    root.bind('<Return>', (lambda event, e=ents: fetch(e)))

    modelLabel = tk.Label(tab1, text='GFET Model')
    modelCombo = ttk.Combobox(tab1, state='readonly')

    modelCombo['values'] = 'Rodriguez', 'Thiele'
    modelCombo.current(0) # First value in list is default
    
    root.bind("<<ComboboxSelected>>")

    modelLabel.grid(row=1, column=0)
    modelCombo.grid(row=1,column=1)
    modelLabel.pack()
    modelCombo.pack()    

    dielecLabel = tk.Label(tab2, text='Relative Permittivity')
    dielecCombo = ttk.Combobox(tab2, state='readonly')
    vals = []
    i = 0
    for i in range(len(dielecs)):
        vals.append(dielecs[i][0] + " (" + str(dielecs[i][1]) + ")")

    dielecCombo['values'] = vals
    dielecCombo.current(0) # First value in list is default
    
    root.bind("<<ComboboxSelected>>")

    dielecLabel.grid(row=1, column=0)
    dielecCombo.grid(row=1,column=1)
    dielecLabel.pack()
    dielecCombo.pack()

    ents2 = makeform(tab2, fields2, params2)
    root.bind('<Return>', (lambda event, e=ents2: fetch(e)))
    
    fig = plt.Figure(figsize=(5,4), dpi=100)
    ax = fig.add_subplot(111)
    ax.set_title('GFET Transfer Characteristic')
    ax.set_ylabel(r'$I_{DS}$ (A)')
    ax.set_xlabel(r'$V_{GS}$ (V)')
    ax.set_aspect(1./ax.get_data_ratio())
    ax.ticklabel_format(axis='both',style='sci', scilimits=(0,0))
      
    b1 = tk.Button(tab1, text='Simulate',
                  command=(lambda e=ents, e2=ents2: calculateIds(e, e2, modelCombo, dielecCombo, ax, canvas)))
    b1.pack(side=tk.LEFT, padx=5, pady=5)
    b2 = tk.Button(tab1, text='Load Sweep', command=loadSweep)
    b2.pack(side=tk.LEFT, padx=5, pady=5)
    b3 = tk.Button(tab1, text='Export CSV', command=(lambda : exportData(V, I)))
    b3.pack(side=tk.LEFT, padx=5, pady=5)
    
    canvas = FigureCanvasTkAgg(fig, master=right)
    canvas.get_tk_widget().pack()
    canvas.draw()

    root.mainloop()
    
