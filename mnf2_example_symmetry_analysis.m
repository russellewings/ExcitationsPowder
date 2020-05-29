% Example scripts to load ISIS neutron training course data on MnF2 using
% Horace, and to simulate and fit using SpinW.
% 
% Assumes that you have already installed Horace and SpinW!
%
% ============
% Russell Ewings - 20/5/2020
%

%% How to examine bonds to determine symmetry-allowed interactions
%
% Especially useful for establishing whether DM is allowed by symmetry, and
% then perhaps whether off-diagonal terms are allowed in the exchange
% tensor


mnf2 = spinw;
mnf2.genlattice('lat_const', [4.87 4.87 3.31], 'angle', [90 90 90]*pi/180, 'sym', 'P 42/m n m');
mnf2.addatom('r', [0 0 0], 'S', S, 'label', 'MMn2', 'color', 'b')
mnf2.gencoupling('maxDistance', 5)

%assign a dummy (zero) exchange to the nth bond (3rd in this example)
mnf2.addmatrix('label','dummy','value',0);
mnf2.addcoupling('mat','dummy','bond',3)

%Now run the getmatrix command to find out about what form of exchange and
%DM coupling is allowed by symmetry. N.B. this can include the existence or
%otherwise of off-diagonal terms in the exchange tensor, in addition to DM.
mnf2.getmatrix('mat','dummy')

%Information about what is allowed by symmetry is printed in the main
%Matlab window

%===

%In th case of MnF2 we find that bonds 1 and 3 do not permit a DM interaction,
%and bond 2 only permits a DM vector of the form [DM1,DM1,0], i.e. no
%z-component and the x and y components equal.
%
%Information is also given about symmetry-allowed differences in the
%diagonal and off-diagonal terms in the exchange tensor
%

%====

%Much more information about symmetry analysis is presented in tutorial no.
%32 on the SpinW website - highly recommended!

%https://spinw.org/tutorials/32tutorial

