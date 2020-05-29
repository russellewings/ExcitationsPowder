function dq=qres(ei,eps,guide)
% FWHH Q resolution in small angle and small fractional energy transfer limit for MERLIN
%
% T.G.Perring 13 August 2008
%
% Use Angular divergence  of incident beam from Monte Carlo simulations by Rob Bewley
%
% dq=qres(ei,eps,guide)


% Widths of components to Q resolution on MERLIN
x1=12;      % moderator-sample distance (m)
x2=2.5;     % sample-detector distance
ws=0.042;   % diameter of sample - treat as cylindical shell
wd=0.025;   % diameter of detector - treat as hlafway between cylinder and cylindrical shell

% Effective full width of moderator (m) - treat as hat function
if ~exist('guide','var')||guide
    % Angular divergence: digitise points on plot from Rob Bewley
    xmc=log([9,30,70,150,400]);%list of energies at which divergence of beam was calculated by Rob
    divmc=[2.0035    1.1849    0.8395    0.6349    0.5070];%list of angular divergence FWHM in degrees at thos energies
    div=interp1(xmc,divmc,log(ei),'spline','extrap');
    wm=div*(pi/180)*x1; % effective moderator full width
else
    wm=0.1;
end

ki=sqrt(ei/2.07214);
kf=sqrt((ei-eps)/2.07214);

% treat each contribution as a hat function - approx resolution

sig_q=ki*(wm/sqrt(12))/x1;
% disp(['         moderator ',num2str(sig_q)])
sig_q=ki*(ws/sqrt(8))/x1;
% disp(['    primary sample ',num2str(sig_q)])
sig_q=kf*(ws/sqrt(8))/x2;
% disp(['  secondary sample ',num2str(sig_q)])
sig_q=kf*(wd/sqrt(12))/x2;
% disp(['          detector ',num2str(sig_q)])
% disp(' ')

var_q=(((wm^2/12)+(ws^2/8))/(x1^2))*ki^2 + (((wd^2/12)+(ws^2/8))/(x2^2)).*kf.^2;
% disp(['           overall ',num2str(sqrt(var_q))])

% Convert variance to FWHH - treat as Gaussian for timebeing
dq=sqrt(var_q)*sqrt(log(256));
% disp(['      overall FWHH ',num2str(dq)])
