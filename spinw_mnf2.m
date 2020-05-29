% Example scripts to load ISIS neutron training course data on MnF2 using
% Horace, and to simulate and fit using SpinW.
% 
% Assumes that you have already installed Horace and SpinW!
%
% ============
% Russell Ewings - 20/5/2020
%

%% Example of using SpinW to plot dispersion and powder average S(Q,E)

% Exchange parameters from: 
% G. G. Low, et al J. Appl. Phys. 35, 998 (1964)
S = 5/2;
J1 = -0.32 * S / 11.6;
J2 = 1.76 * S / 11.6;
D = 0;%NB multiplied by zero, so no anisotropy

J1=J1/1.2; J2=J2/1.2; D=D/1.2;

% Setup spinw model
mnf2 = spinw;
mnf2.genlattice('lat_const', [4.87 4.87 3.31], 'angle', [90 90 90]*pi/180, 'sym', 'P 42/m n m');
mnf2.addatom('r', [0 0 0], 'S', S, 'label', 'MMn2', 'color', 'b')
mnf2.gencoupling('maxDistance', 5)
mnf2.addmatrix('label', 'J1', 'value', J1, 'color', 'red');
mnf2.addmatrix('label', 'J2', 'value', J2, 'color', 'green');
mnf2.addcoupling('mat', 'J1', 'bond', 1)
mnf2.addcoupling('mat', 'J2', 'bond', 2)
mnf2.addmatrix('label', 'D', 'value', diag([0 0 D]), 'color', 'black');
mnf2.addaniso('D')
mnf2.genmagstr('mode', 'direct', 'S', [0 0; 0 0; 1 -1])

% Plots magnetic structure and bonds
plot(mnf2)

% Calculates the powder spectrum
spectrum = mnf2.powspec([0:0.01:4.5], 'Evect', [0:0.05:12], 'nRand', 1e3, 'fibo', false, 'hermit', true,...
    'formfact',true);
figure;
sw_plotspec(spectrum,'mode','color')
caxis([0 0.5])
colormap parula


% Convolutes with a fixed linewidth and applies kinematic limits
spectrum = sw_instrument(spectrum, 'dE', 0.3, 'dQ', 0.01, 'ThetaMin', 2.5, 'Ei', 12,'formfact',true);
figure;
sw_plotspec(spectrum,'mode','color')
caxis([0 0.3])
colormap parula


%==============
%Plot crystal spin wave dispersion
Qcorner = {[1/2 1/2 0] [0 0 1] [0 0 0] [1/2 1/2 0] [1/2 1/2 1/2] [1/2 0 1/2] [0 0 0] [1/2 0 0] 200};
spec = mnf2.spinwave(Qcorner,'hermit',false);

figure
sw_plotspec(spec,'mode','disp','imag',true,'colormap',[0 0 0],'colorbar',true)

