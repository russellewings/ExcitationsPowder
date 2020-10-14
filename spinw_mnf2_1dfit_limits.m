function yout=spinw_mnf2_1dfit_limits(E,pars,parlims,Qrange,ok,Ei,dE,dQ,s)

% Fitting a 1d cut along the energy axis.
% Enforce limits one or more parameters using the additional input
% parameter "parlims" - this is a cell array with the same number of
% elements as the pars array. Each element of the cell array should be
% either an empty array, given by [], or a 2 element array of the form
% [lo_limit,hi_limit]. The nth element of the cell array corresponds to the
% nth element of the parameter array

%Perform some checks:
if numel(pars)~=numel(parlims)
    error('Check number of fit parameters and number of limits are consistent');
end
if ~iscell(parlims)
    error('parlims input must be a cell array');
end
for i=1:numel(parlims)
    if isempty(parlims{i}) || isnumeric(parlims{i})
        %all good
    else
        error('Check that all elements of parlims are either empty arrays, or are 2 element numeric arrays');
    end
    if isnumeric(parlims{i}) && ~isempty(parlims{i}) && numel(parlims{i})~=2
        error('Check that all elements of parlims are either empty arrays, or are 2 element numeric arrays');
    end
end
%If we get to this point the input format is probably correct...

%Now implement setting the parameters to be within the limits. This does
%then mean that the output of the fitting programme needs to be checked -
%see main code that calls this for explanation.

S=5/2;
for i=1:numel(pars)
    if ~isempty(parlims{i})
        p(i)=parlims{i}(1) + (0.5.*(sin(pars(i)) + 1)).*(parlims{i}(2) - parlims{i}(1));
    else
        p(i)=pars(i);
    end
end
    
scalefac = p(1);
J=p(2:3);
D=p(4);
bg=p(5);

%Ensure reasonable number of Q bins between limits for the cut. May wish to
%hand-craft this. Generally want to make QQ have a step size comparable to
%the 2d plot of the signal vs |Q| and E.
QQ=linspace(Qrange(1),Qrange(2),30);

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
    mnf2powspec=mnf2.powspec_ran(QQ','Evect',unique(E),...
        'binType','cbin','nRand',1000,'hermit',true,...
        'formfact',true,'s_rng',s);%note that you can change the number of random Q points
    %A smaller number gives a noisier output, but faster evaluation. For
    %this 1d cut can afford to make it larger than for the 2d case.

    %Give a file containing Nx2 matrix giving Etrans and dE. From Pychop for Ei=11meV 240/120Hz
    mnf2powspec = sw_instrument(mnf2powspec,'dE',dE,...
        'Ei',Ei,'dQ',dQ);

    yout=abs(scalefac).*mnf2powspec.swConv';
    
    %Extra catch for very small imaginary values of S(Q,w):
    rr=real(yout);
    ii=imag(yout);
    mr=max(max(abs(rr)));
    mi=max(max(abs(ii)));
    if mi<1e-5*mr
        yout=rr;
    else
        disp('Warning: Significant imaginary values for S(Q,w) in SpinW calculation');
        error('Significant imaginary values for S(Q,w) in SpinW calculation');
    end
    
    %Do the integration:
    yout=sum(yout,1);
    
    yout=yout+bg;
    
    yout=yout(ok);
    
    %Can be odd cases when small number of additional points come from sim
    %as NaN. In this case replace them with bg:
    f=isnan(yout);
    yout(f)=bg;
    
catch
    %deal with the case that spinwave calc failed - e.g. exchanges
    %inconsistent with known magnetic structure. Returns a signal array
    %that is so massive that it will ensure the fitting algorithm avoids
    %such places
    yout=ones(1,numel(unique(E))).*1e12;
    yout=yout(ok);
end