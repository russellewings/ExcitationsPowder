function yout=spinw_mnf2_1dfit_pso(E,pars)

% Fitting a 1d cut along the energy axis, set up for particle swarm
% optimization. Because of limitations of the SpinW particle optimizer we
% are only allowed to pass in the "x" value and the input parameters. So we
% have to hard code in the other parts from similarly named functions. We
% can fix parameters in our fitting with the pso optimizer by setting the
% upper and lower bounds on their values to be identical.
%

% Options we have to hard code and/or keep fixed:
%Qrange,ok,Ei,dE,dQ,s

S=5/2;

%"Real" fitting parameters
scalefac = pars(1);
J=pars(2:3);
D=pars(4);
bg=pars(5);

%Parameters we need that we should never fit
Q1=pars(6); Q2=pars(7);%these define the Q range
Ei=pars(8);
dE=pars(9);
dQ=pars(10);
rnseed=pars(11);%integer seed for the random number generator, so that the
%same set of random numbers are used in the evaluation of the powder
%average every time.


%Ensure reasonable number of Q bins between limits for the cut. May wish to
%hand-craft this. Generally want to make QQ have a step size comparable to
%the 2d plot of the signal vs |Q| and E.
QQ=linspace(Q1,Q2,30);

%Random number setup:
s=rng(rnseed);

% Setup spinw model
mnf2 = spinw;
mnf2.genlattice('lat_const', [4.87 4.87 3.31], 'angle', [90 90 90]*pi/180, 'sym', 'P 42/m n m');
mnf2.addatom('r', [0 0 0], 'S', S, 'label', 'MMn2', 'color', 'b')
mnf2.gencoupling('maxDistance', 5)
mnf2.addmatrix('label', 'J1', 'value', J(1), 'color', 'red');
mnf2.addmatrix('label', 'J2', 'value', J(2), 'color', 'green');
mnf2.addcoupling('mat', 'J1', 'bond', 1)
mnf2.addcoupling('mat', 'J2', 'bond', 2)
mnf2.addmatrix('label', 'D', 'value', diag([0 0 D]), 'color', 'black');
mnf2.addaniso('D')
mnf2.genmagstr('mode', 'direct', 'S', [0 0; 0 0; 1 -1])

try
    % Powder average spin waves:
    mnf2powspec=mnf2.powspec_ran(QQ','Evect',unique(E)',...
        'binType','cbin','nRand',1000,'hermit',true,...
        'formfact',true,'s_rng',s);%note that you can change the number of random Q points
    %A smaller number gives a noisier output, but faster evaluation. For
    %this 1d cut can afford to make it larger than for the 2d case.

    %Give a file containing Nx2 matrix giving Etrans and dE. From Pychop for Ei=11meV 240/120Hz
    mnf2powspec = sw_instrument(mnf2powspec,'dE',dE,...
        'Ei',Ei,'dQ',dQ);

    yout=abs(scalefac).*mnf2powspec.swConv';
    %Do the integration:
    yout=sum(yout,1);
    
    yout=yout+bg;
    
    %We skip this bit for particle swarm optimisation
    %yout=yout(ok);
    
    %Can be odd cases when small number of additional points come from sim
    %as NaN. In this case replace them with bg:
    f=isnan(yout);
    yout(f)=bg;
    
    %Convert output back to correct size for pso:
    yout=yout';
    
catch
    %deal with the case that spinwave calc failed - e.g. exchanges
    %inconsistent with known magnetic structure. Returns a signal array
    %that is so massive that it will ensure the fitting algorithm avoids
    %such places
    yout=ones(1,numel(unique(E))).*1e12;
    %yout=yout(ok);
end