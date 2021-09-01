function [Tinf,pinf,rhoinf,cpinf,viscinf,R,cvinf,gamma,a_sound,qinf,Mach,Re,Kn,Cpmax]=AtmosphereCuriosity(h,Vinf,Diameter)
load MCDCuriosityAtmosphere.mat
if h>250000
    Tinf=MCDCuriosityAtmosphere{end,2};     % Freestream temperature [K]
    pinf=MCDCuriosityAtmosphere{end,3};     % Freestream pressure [k]
    rhoinf=MCDCuriosityAtmosphere{end,4};   % Freestream density [kg/m^3]
    cpinf=MCDCuriosityAtmosphere{end,5};    % Specific heat at constant pressure
    viscinf=MCDCuriosityAtmosphere{end,6};  % Dynamic viscosity [Ns/m2]
else
    Tinf=interp1(MCDCuriosityAtmosphere{:,1},MCDCuriosityAtmosphere{:,2},h,'spline');
    pinf=interp1(MCDCuriosityAtmosphere{:,1},MCDCuriosityAtmosphere{:,3},h,'spline');
    rhoinf=interp1(MCDCuriosityAtmosphere{:,1},MCDCuriosityAtmosphere{:,4},h,'spline');
    cpinf=interp1(MCDCuriosityAtmosphere{:,1},MCDCuriosityAtmosphere{:,5},h,'spline');
    viscinf=interp1(MCDCuriosityAtmosphere{:,1},MCDCuriosityAtmosphere{:,6},h,'spline');
end
% freestream properties
R=pinf/(rhoinf*Tinf);               % Gas constant R
cvinf=cpinf-R;                      % Specific heat at constant volume Cv
gamma=cpinf/cvinf;                  % Ratio of specific heats
a_sound=sqrt(gamma*R*Tinf);         % Speed of sound
qinf=0.5*rhoinf*Vinf^2;             % Dynamic pressure
Mach=Vinf/a_sound;                  % Mach number
Re=rhoinf*Vinf*Diameter/viscinf;    % Reynolds Number
Kn=Mach/Re*sqrt(pi*gamma/2);        % Knudsen Number
Cpmax=(2/(Mach^2*gamma))*((((gamma+1)^2*Mach^2)/(4*gamma*Mach^2-2*(gamma-1)))^(gamma/(gamma-1))*((1-gamma+2*gamma*Mach^2)/(gamma+1))-1);
return
