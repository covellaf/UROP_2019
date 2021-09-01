function vardot=Trajectory(t,var)
global x y z TrajResults 
p = readtable('GeneralData.csv.xlsx', 'ReadVariableNames',false);  
rMars = p.Var2(1);
omegaMars = p.Var2(2);
muMars = p.Var2(3);
kMars = p.Var2(4);
J2Mars = p.Var2(5);
CMars = p.Var2(6);
aMars = p.Var2(7);
bMars = p.Var2(8);
Diameter = p.Var2(44);
RNose = p.Var2(45);
mass = p.Var2(48);
RNoseCone = p.Var2(49);
RRib = p.Var2(50);
NRibs = p.Var2(51);
RibAngle = p.Var2(52);
RefArea = p.Var2(54);
RPayload = p.Var2(55);
HeightPayload = p.Var2(56);
HSthickness = p.Var2(57);
RibWidth = p.Var2(58);
RibDepth = p.Var2(59);
RibWallThickness = p.Var2(60);
EndRibWidth = p.Var2(62);
EndRibDepth = p.Var2(63);
EndRibWallThickness = p.Var2(64);
SupportStrutLoc = p.Var2(68);
RibLength = p.Var2(70);
E = p.Var2(71);
HSdensity = p.Var2(74);
SphereConeAngleRad = p.Var2(79);                  
NoseHeight = p.Var2(80);                   
NoseLength = p.Var2(81);
ConeVertex = p.Var2(87);
RibAngleRad = p.Var2(88);
nopNose = p.Var2(94);
nopFrustum = p.Var2(95);
nopRibLength = p.Var2(96);
nopRibWidth = p.Var2(97);
%-----------------------------------------------------------------------------
r=var(1); lat=var(2); long=var(3);
u=var(4); v=var(5); w=var(6);
e0=var(7); e1=var(8); e2=var(9); e3=var(10);
omegaX=var(11); omegaY=var(12); omegaZ=var(13);
% Calculate other trajectory variables
h=r-rMars; 
Vinfgc=[u;v;w];
Vinf=sqrt(u^2+v^2+w^2);
FPA=asin(w/Vinf);
azi=atan2(v,u);
roll=atan2(2*e0*e1+2*e2*e3,e0^2-e1^2-e2^2+e3^2);
pitch=-asin(2*(e1*e3-e0*e2));                       % 1,2,3 Quaternion to Euler angle conversion (Diebel, 2006)
yaw=atan2(2*e1*e2+2*e0*e3,e0^2+e1^2-e2^2-e3^2);
omega=[omegaX; omegaY; omegaZ];
% Set up rotation matrices
Ggd2b=[e0^2+e1^2-e2^2-e3^2, 2*(e1*e2+e0*e3), 2*(e1*e3-e0*e2);2*(e1*e2-e0*e3), e0^2-e1^2+e2^2-e3^2, 2*(e0*e1+e2*e3);...
    2*(e1*e3+e0*e2), 2*(e2*e3-e0*e1), e0^2-e1^2-e2^2+e3^2];         % Transformation matrix from geodetic to body frame (Karlgaard, 2006)
Ggc2gd=[1 0 0;0 1 0;0 0 1];     % Assuming Mars is spherical, there is no difference between geocentric and geodetic
Ggc2b=Ggc2gd*Ggd2b;
Vinfb=Ggc2b*Vinfgc;
vx=Vinfb(1); vy=Vinfb(2); vz=Vinfb(3);
beta=atan2(vy,(sqrt(vx^2+vz^2)));
alpha=atan2(vz,vx);
bank=atan2(2*(e0*e1+e2*e3)+sin(beta)*sin(-FPA),((e0^2-e1^2+e2^2-e3^2)*cos(azi)-(2*(e1*e2-e0*e3))*sin(azi))*cos(-FPA));    % From Wagner (1970) eqn. 229
Vinfwind=[Vinf*cos(beta)*cos(alpha); Vinf*sin(beta); Vinf*cos(beta)*sin(alpha)];
alpha_total=acos(cos(alpha)*cos(beta));     % Total angle of attack [rad] - used in Airbus reference frame
phi_a=atan2(tan(beta),sin(alpha));          % Aerodynamic roll angle [rad] - used in Airbus reference frame
% Define atmospheric properties--------------------------------------------
[Tinf,pinf,rhoinf,cpinf,viscinf,R,cvinf,gamma,a_sound,qinf,Mach,Re,Kn,Cpmax]=AtmosphereCuriosity(h,Vinf,Diameter);
% Reset geometry to initial state for next bending
[x,y,z]=Geometry2(RNoseCone,RNose,RRib,NRibs,SphereConeAngleRad,RibAngleRad... line2
    ,nopNose,nopFrustum,nopRibLength,nopRibWidth,ConeVertex,NoseHeight,NoseLength);
% Calculate aerodynamic coefficients and angular accelerations
[Cl,Cd,LDRatio,Croll,Cpitch,Cyaw,AngAccels]=Aerodynamics(x,y,z,Vinfwind,Cpmax,qinf,pinf,alpha,omega... line 2
    ,mass,HSdensity,HSthickness,Diameter,RPayload,HeightPayload,RefArea,RibAngle,E,RibWidth,RibDepth,RibWallThickness... line 3
    ,EndRibWidth,EndRibDepth,EndRibWallThickness,SupportStrutLoc,RibLength,NRibs,nopNose,nopFrustum,nopRibLength,nopRibWidth);
%--------------------------------------------------------------------------
omegaXdot=AngAccels(1);
omegaYdot=AngAccels(2);
omegaZdot=AngAccels(3);
% Calculate heating values
qconv=kMars*sqrt(rhoinf/RNose)*Vinf^3/10000;  % Sutton-Graves relation for convective heat flux at stagnation point [W/cm2]
if Vinf>6500
    load TauberSuttonVelocityFunction
    TSv=interp1(TauberSuttonVelocityFunction{:,1},TauberSuttonVelocityFunction{:,2},Vinf,'spline');
    qrad=CMars*RNose^aMars*rhoinf^bMars*TSv;    % Tauber-Sutton relation for radiative heat flux at stagnation point [W/cm2]
else
    qrad=0;
end
% Calculate aerodynamic forces and decelerations
D=Cd*qinf*RefArea;                  % Drag force [N]
L=Cl*qinf*RefArea;                  % Lift force [N]
BC=mass/(Cd*RefArea);
Du=-D*cos(FPA)*cos(azi);    % Component of drag force in N direction [N]
Dv=-D*cos(FPA)*sin(azi);    % Component of drag force in E direction [N]
Dw=-D*sin(FPA);             % Component of drag force in D direction [N]
Lu=L*sin(FPA)*cos(azi);     % Component of lift force in N direction [N]
Lv=L*sin(FPA)*sin(azi);     % Component of lift force in E direction [N]
Lw=-L*cos(FPA);             % Component of lift force in D direction [N]
au=(Du+Lu)/mass;            % Component of aerodynamic acceleration in N direction [m/s2]
av=(Dv+Lv)/mass;            % Component of aerodynamic acceleration in E direction [m/s2]
aw=(Dw+Lw)/mass;            % Component of aerodynamic acceleration in D direction [m/s2]
a=sqrt(au^2+av^2+aw^2);     % Total aerodynamic acceleration
gload=a/9.81;                                     % Total g-load imparted to vehicle
J=3/2*J2Mars*rMars^2;                             % Wagner (1970) definition
gu=-(muMars*J*sin(2*lat)/(r^4));                  % Component of gravitational acceleration in N direction [m/s2] Wagner (1970)
gv=0;                                             % Component of gravitational acceleration in E direction [m/s2]
gw=(muMars/r^2)-(muMars*(2-3*cos(lat).^2)/(r^4)); % Component of gravitational acceleration in D direction [m/s2]
% Run equations of motion
rdot=-w;                    % Rate of change of r [m/s] (Karlgaard, 2009)
latdot=u/r;                 % Rate of change of latitude [rad/s]
longdot=v/(r*cos(lat));     % Rate of change of longitude [rad/s] POSITIVE FOR RETROGRADE
udot=au+(1/r)*(u*w-v^2*tan(lat))+gu-2*omegaMars*v*sin(lat);%-omegaMars^2*r*sin(2*lat)/2;    % Acceleration in N direction [m/s2] (Karlgaard, 2009)
vdot=av+(1/r)*(u*v*tan(lat)+v*w)+gv+2*omegaMars*(w*cos(lat)+u*sin(lat));    % Acceleration in E direction [m/s2]
wdot=aw-(1/r)*(u^2+v^2)+gw-2*omegaMars*v*cos(lat);%-omegaMars^2*r/2*(1+cos(2*lat));             % Acceleration in D direction [m/s2]
e0dot=0.5*[-e1 -e2 -e3]*([omegaX; omegaY; omegaZ]-(1/r)*Ggc2b*[v; -u; -v*tan(lat)]-Ggc2b*[omegaMars*cos(lat); 0; -omegaMars*sin(lat)]);
e1dot=0.5*[e0 -e3 e2]*([omegaX; omegaY; omegaZ]-(1/r)*Ggc2b*[v; -u; -v*tan(lat)]-Ggc2b*[omegaMars*cos(lat); 0; -omegaMars*sin(lat)]);
e2dot=0.5*[e3 e0 -e1]*([omegaX; omegaY; omegaZ]-(1/r)*Ggc2b*[v; -u; -v*tan(lat)]-Ggc2b*[omegaMars*cos(lat); 0; -omegaMars*sin(lat)]);
e3dot=0.5*[-e2 e1 e0]*([omegaX; omegaY; omegaZ]-(1/r)*Ggc2b*[v; -u; -v*tan(lat)]-Ggc2b*[omegaMars*cos(lat); 0; -omegaMars*sin(lat)]);
vardot=[rdot; latdot; longdot; udot; vdot; wdot; e0dot; e1dot; e2dot; e3dot; omegaXdot; omegaYdot; omegaZdot];
TrajResults=[TrajResults;t,h,Vinf,Mach,rad2deg(FPA),rad2deg(azi),rad2deg(lat),rad2deg(long),u,v,w,BC,gload,rad2deg(alpha),rad2deg(beta),...
    e0,e1,e2,e3,rad2deg(pitch),rad2deg(yaw),rad2deg(roll),rad2deg(bank),rad2deg(alpha_total),rad2deg(phi_a),rad2deg(omegaX),rad2deg(omegaY),...
    rad2deg(omegaZ),omegaXdot,omegaYdot,omegaZdot,Tinf,pinf,rhoinf,cpinf,viscinf,R,cvinf,gamma,a_sound,qinf,qconv,qrad];
return