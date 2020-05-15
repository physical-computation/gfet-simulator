# gfet-simulator
Python Program to Model GFET Characteristics

Author: Nathaniel Tye (njt48@cam.ac.uk)
Date: 01.05.2020

Purpose:
This program is designed as a tool allowing people to play around with different parameters of GFETs,
e.g. channel dimensions, oxide thickness, operating temperature, and to test these devices under a range
of bias conditions. Currently, Rodriguez et al. [1] and Thiele et al.'s [2] models have been implemented.

Motivation:
GFETs are being used in the research of several people in the group, and potentially by new students as well,
but they aren't off the shelf components, and can take hours to make. As they are also very sensitive devices,
a way for people to be able to experiment with GFETs, without the need for fabrication or the risk of damaging
devices is useful. There are currently very few publicly-available tools out there that allow this and many 
models in literature, from my reading and at the time of writing, seem to exist either as closed-source, 
limited-access (and limited-features) applications, or just as equations or proposed SPICE implementations.
The intention here is then to save people having to spend time trying to understand the derivations and 
differences in existing GFET models, and instead be able to delve straight in and get on with their research.


Current Features:
- User can specify physical device parameters (channel dimensions, gate oxide thickness,
  Carrier density, dielectric material, carrier mobility)
- Simulates both I-V and Transfer characteristics with user-defined sweeps
- Can do linear, dual-linear and logarithmic sweeps of gate-source and drain-source
  Voltages
- Can export CSV files of the measured transfer and IV characteristics
- Can load in a user-defined sweep from excel (limited functionality at present)

Features to add:
- Additional models (e.g. high-frequency, multi-gate [e.g. dual top gate, back gate] )
- Include model parameters and other information in exported CSV files
- Improve user-defined sweep loading
- Plot additional characteristics (e.g. conductivity, capacitance, frequency response...)
    -- Include parameters such as temperature to this end
- Possibly even export a SPICE model of the device tested, helping 
  facilitate rapid circuit-level simulation

Other general program functionality todo:
- Input validation & general stability
- Aesthetics

References:
[1] S. Rodriguez et al., "A Comprehensive Graphene FET Model for Circuit Design," 
    in IEEE Transactions on Electron Devices, vol. 61, no. 4, pp. 1199-1206, April 2014. (doi.org:/10.1109/TED.2014.2302372)

[2]
