% Example scripts to load ISIS neutron training course data on MnF2 using
% Horace, and to simulate and fit using SpinW.
% 
% Assumes that you have already installed Horace and SpinW!
%
% ============
% Russell Ewings - 20/5/2020
% (and Richard Dixey - 30/3/2021)

%% Simulation with Spinw

%First step is to convert the data (from mnf2_example_plotting.m) into yet another format!
mnf2_IX=IX_dataset_2d(mnf2_cut.p{1},mnf2_cut.p{2},mnf2_cut.s,sqrt(mnf2_cut.e));

%Notice that this has zero intensity whereas before there were NaNs to
%specify parts of QE space where there are no detectors. Change this back
%so we know where we measured.

mnf2_IX.signal(mnf2_IX.error==0 & mnf2_IX.signal==0)=NaN;
plot(mnf2_IX)
%shold look the same as before

%Use the Horace multifit tool to simulate and/or fit a spinwave mode. See
%spinw_mnf2.m for details of doing a regular simulation, spinw_mnf2_fit.m
%for use with these routines

%First just evaluate the cross section using the known correct values for
%J:
scalefac=1;%intensity scale factor
J=[-0.0575,0.3161];%exchange (correct values used here)
D=0;%single ion anisotropy
bg=10;%background

%Tell the function which bits of QE space are OK (i.e. are not NaN):
%Includes an important extra consideration, spotted by Richard Dixey of
%QMUL, of what happens if an entire row or column of the signal matrix is
%NaN. To get around this bug the first 3 lines are necessary:
nonZeroRows = find(all(isnan(mnf2_IX.signal),2)); nonZeroCols = find(all(isnan(mnf2_IX.signal),1));
mnf2_IX.signal(nonZeroRows,1)=1e-6; mnf2_IX.error(nonZeroRows,1)=1e6;
mnf2_IX.signal(1,nonZeroCols)=1e-6; mnf2_IX.error(1,nonZeroCols)=1e6;
%This replaces one element with a very small signal with a very big error,
%just to avoid any weight being given to these points in fits.

%Now find which points are "OK" or not
ok=~isnan(mnf2_IX.signal);

%Instrument settings
Ei=12;
dE=0.36;%Check with your instrument scientist! This can also be the name of
%a  file that gives the resolution as a function of energy transfer
dQ=0.04;
s=rng(1);%random number seed, so that you get the same result each time
%(the SpinW powder averaging routine takes the average of many Q points
%with the same value of |Q|, randomly chosen, for the powder average - by
%seeding the random number generator we ensure that these Q points are the
%same every time we run the simulation, which aids fitting)

%===========================
% <REALLY IMPORTANT POINT>:
% In the ExcitationsPowder routines that you downloaded, you must copy the file "powspec_ran.m"
% to the folder in your SpinW installation .../swfiles/@spinw (this is where the regular powspec.m
% file lives - check by typing in Matlab "which powspec"). The routine ensures that the same random 
% number seed is used every time - see above comment for details
%===========================

[wfit,fitdata]=multifit(mnf2_IX,@spinw_mnf2_fit,{[scalefac,J,D,bg],ok,Ei,dE,dQ,s},...
    [1,1,1,0,1],'evaluate');
%
% We use the 'evaluate' keyword so that a fit is not done, simply the
% function is evaluated using our initial guess. If we omitted this keyword
% a fit would be done, with the first 3 parameters (scale factor and Js)
% being fitted, but D fixed and background fitted: [1,1,1,0,1] - 1 for free, 0 for fixed. 

%compare data and simulation:
plot(mnf2_IX); ly 0 12; lz 0 80; keep_figure;
plot(wfit); ly 0 12; lz 10 10.6; keep_figure

%We can see that we can get an excellent agreement if we use the right
%parameters!