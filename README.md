# GFET Lab
## A Python Program to Model GFET Characteristics

Author: Nathaniel Tye (njt48@cam.ac.uk)
Date: 01.05.2020
Updated: 28.06.2022

## Purpose:
This program is designed as a tool allowing people to experiment with different graphene field-effect transistor (GFET) parameters,
e.g. channel dimensions, oxide thickness, operating temperature, and to test these devices under a range
of bias conditions. We have implemented Jimenez's model [1], but we encourage community integration of other models.

Please see Usage.md for a detailed guide on use of the software.

The accompanying preprint paper is available [here](https://doi.org/10.48550/arxiv.2206.13239).

To cite this article and tool, use the following:

>@misc{https://doi.org/10.48550/arxiv.2206.13239,
>  doi = {10.48550/ARXIV.2206.13239},
>  url = {https://arxiv.org/abs/2206.13239},
> author = {Tye, Nathaniel J. and Tadbier, Abdul Wadood and Hofmann, Stephan and Stanley-Marbell, Phillip},
>  keywords = {Mesoscale and Nanoscale Physics (cond-mat.mes-hall), Applied Physics (physics.app-ph), FOS: Physical sciences, FOS: Physical sciences},
>  title = {GFET Lab: A Graphene Field-Effect Transistor TCAD Tool},
>  publisher = {arXiv},
>  year = {2022},
>  copyright = {arXiv.org perpetual, non-exclusive license}
>}


## Dependencies:
This tool is written in Python 3.9 and depends on several external libraries for correct operation. Please make sure you have the following installed before use:

* NumPy (https://numpy.org/)
* SciPy (https://www.scipy.org/)
* MatPlotLib (https://matplotlib.org/stable/index.html)
* Tkinter (Included with Python 3.7+)

Prebuilt executables for Mac and Windows are available in the respective folders. Note that the Windows version requires "Dielectrics.txt" to be located in the same directory.

## Motivation:
GFETs are being used in the research of several people in the group, and potentially by new students as well,
but they aren't off the shelf components, and can take hours to make. As they are also very sensitive devices,
a way for people to be able to experiment with GFETs, without the need for fabrication or the risk of damaging
devices is useful. There are currently very few publicly-available tools out there that allow this and many 
models in literature, from my reading and at the time of writing, seem to exist either as closed-source, 
limited-access (and limited-features) applications, or just as equations or proposed SPICE implementations.
The intention here is then to save people having to spend time trying to understand the derivations and 
differences in existing GFET models, and instead be able to delve straight in and get on with their research.

## Current Features:
* User can specify physical device parameters (channel dimensions, gate oxide thickness,
  Carrier density, dielectric material, carrier mobility)
* Simulates both I-V and Transfer characteristics with user-defined sweeps
* Can do linear, dual-linear and logarithmic sweeps of gate-source and drain-source
  voltages
* Can export CSV files of the measured transfer and IV characteristics
* Can load in a user-defined sweep from a CSV file for both gate and drain voltages 
  (Currently defaults back-gate voltages to 0 V for external sweeps.)
* Can export SPICE models

## Features to Add:
* Additional models (e.g. multi-gate), alternative materials (e.g. MoS2)...
* Include model parameters and other information in exported CSV files
* Plot additional characteristics (e.g. conductivity, capacitance...)
* 'Help' tab/window
* Schematic of device in current model

## Other Functionality to Implement:
* Ensure stability
* Aesthetics ('simulation running' label at the bottom, or some sort of popup window that closes when finished)

## Known Issues:
* Transconducatance and frequency calculations are placeholders and need to be implemented.
* Mukherjee's model [2] does not currently work, but is instead a placeholder.
* Release v0.1-2d-conf: Exported SPICE models can only consider a voltage input to the drain terminal, the next version will correct this, but it is not an issue for basic I-V and transfer characteristic sweeps.
     
## References:
[1] D. Jimenez, "Explicit Drain Current, Charge and Capacitance Model of Graphene Field-Effect Transistors," in IEEE Transactions on Electron Devices, vol. 58, no. 12, pp. 4377-4383, Dec. 2011, (https://doi.org/10.1109/TED.2011.2168960)

[2] C. Mukherjee, J. Aguirre-Morales, S. Frégonèse, T. Zimmer and C. Maneux, "Versatile Compact Model for Graphene FET Targeting Reliability-Aware Circuit Design," in IEEE Transactions on Electron Devices, vol. 62, no. 3, pp. 757-763, March 2015, doi: 10.1109/TED.2015.2395134.
