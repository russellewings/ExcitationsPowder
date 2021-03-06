# ExcitationsPowder
Selection of powder data analysis routines, provided without warranty

Description of routines:

Algorithms:

powspec_ran - <IMPORTANT> this must be copied to your spinw directory .../swfiles/@spinw/. The routine ensures that powder spectra are calculated by SpinW using the
same random number seed every time, which is crucial for fitting to work

qres - a utility routine to calculate the transverse q resolution in the small angle approximation. It is written just for MERLIN right now, but changing the value 
x2 in the code (secondary flight path) will make results roughly OK for MAPS or MARI also.

---
Template scripts:

mnf2_example_plotting - demonstrates use of Horace and associated tools for plotting powder data, and taking cuts

mnf2_example_simulating_spinw - demonstrates use of Horace and SpinW when combined for simulating the powder cross-section due to spin waves

mnf2_example_fitting_spinw - as above, with a description of different fitting strategies for use of SpinW with powder data

mnf2_example_symmetry_analysis - an introduction to using SpinW's in-built symmetry analysis routines, e.g. to determine which bonds can support a DM interaction.

---
Template routines for analyis, that are called by the template scripts. Some hand-crafting required for use with your own data:

spinw_mnf2 (simulation of MnF2 data [2d])

spinw_mnf2_fit (fitting of MnF2 data [2d])

spinw_mnf2_1dfit (fitting of MnF2 data [1d cuts])

spinw_mnf2_1dfit_limits (as above with limits applied on fit parameters)

spinw_mnf2_1dfit_pso (fitting of MnF2 data [1d cuts] using particle swarm optimisation)
