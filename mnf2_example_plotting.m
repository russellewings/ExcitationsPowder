% Example scripts to load ISIS neutron training course data on MnF2 using
% Horace, and to simulate and fit using SpinW.
% 
% Assumes that you have already installed Horace and SpinW!
%
% ============
% Russell Ewings - 20/5/2020
%

%% Get data into format that can be plotted and fitted with Horace routines

%Data file
spe_file='C:\Russell\Harry\CaFe2O4\MAR18301_Ei12.00meV.nxspe';

%If in nxspe format should not need "par" file, which specifies detector
%positions, as already contained in the data file. Use blank string here
par_file='';

%Name of your choice for Horace sqw file
sqw_mnf2='C:\Russell\Harry\CaFe2O4\MnF2.sqw';

%Tell Horace we are using a direct geometry spectrometer
emode=1;

%Set incident energy
efix=12.0;

%Make the sqw file that we then work with later
gen_sqw_powder_test (spe_file, par_file, sqw_mnf2, efix, emode);

%=====

%Take a cut that encompasses the entire data range. Notice the non-standard
%(for Horace) method of *not* specifying a projection axis

mnf2_cut=cut_sqw(sqw_mnf2,0.03,0.12,'-nopix');
%In the above, 0.03 specifies the bin width in |Q|, 0.12 specifies the bin
%width in energy transfer. '-nopix' means we don't bother retaining
%detector pixel information (see Horace documentation for further details)

%Plot this:
plot(mnf2_cut);

%Use Horace commands to change colour scale and axes limits (see Horace
%manual for further plotting options)
lz 0 80;%colour scale
lx 0 4.5;%x-axis limits
ly 0 10;%y-axis limits

%Take a 1d cut and plot it:
mnf2_cut2=cut_sqw(sqw_mnf2,[1.2,1.4],[1,0.12,9],'-nopix');
%Here we've specified to integrate the signal between 1.2<Q<1.4, and make a
%cut along the energy axis from 1 to 9 in steps of 0.12meV

plot(mnf2_cut2)

%If all you wish to do is plot the data, this should be enough to get you
%started.