function [Cl,Cd,LDRatio,Croll,Cpitch,Cyaw,AngAccels]=Aerodynamics(x,y,z,Vinfwind,Cpmax,qinf,pinf,alpha,omega... line 2
    ,mass,HSdensity,HSthickness,Diameter,RPayload,HeightPayload,RefArea,RibAngle,E,RibWidth,RibDepth,RibWallThickness... line 3
    ,EndRibWidth,EndRibDepth,EndRibWallThickness,SupportStrutLoc,RibLength,NRibs,nopNose,nopFrustum,nopRibLength,nopRibWidth)
global PsFs XForces YForces ZForces
COMPL=[-(0.27*Diameter); 0; 0]; % Define payload centre of mass location
xvector=[1 0 0];                % Define x-axis vector
yvector=[0 1 0];                % Define y-axis vector
zvector=[0 0 1];                % Define z-axis vector
[Lx,Ly] = size(x); 
%initialising matrices-----------------------------------------------------
Areas=zeros(Lx-1,Ly-1); Volumes=zeros(Lx-1,Ly-1); Masses=zeros(Lx-1,Ly-1);
ThetasX=zeros(Lx-1,Ly-1); ThetasY=zeros(Lx-1,Ly-1); ThetasZ=zeros(Lx-1,Ly-1); ThetasV=zeros(Lx-1,Ly-1);
HSAngle=zeros(Lx-1,Ly-1);
Cps=zeros(Lx-1,Ly-1);
MaxPressures=zeros(Lx-1,Ly-1);
normalx=zeros(Lx-1,Ly-1); normaly=zeros(Lx-1,Ly-1); normalz=zeros(Lx-1,Ly-1);
PressuresXloc=zeros(Lx-1,Ly-1); PressuresYloc=zeros(Lx-1,Ly-1); PressuresZloc=zeros(Lx-1,Ly-1);
t1locx=zeros(Lx-1,Ly-1); t1locy=zeros(Lx-1,Ly-1); t1locz=zeros(Lx-1,Ly-1);
t2locx=zeros(Lx-1,Ly-1); t2locy=zeros(Lx-1,Ly-1); t2locz=zeros(Lx-1,Ly-1);
NormalForces=zeros(Lx-1,Ly-1);
ForcesX=zeros(Lx-1,Ly-1); ForcesY=zeros(Lx-1,Ly-1); ForcesZ=zeros(Lx-1,Ly-1); ForcesV=zeros(Lx-1,Ly-1);
MassMomentsX=zeros(Lx-1,Ly-1); MassMomentsY=zeros(Lx-1,Ly-1); MassMomentsZ=zeros(Lx-1,Ly-1);
ForceMomentsX=zeros(Lx-1,Ly-1); ForceMomentsY=zeros(Lx-1,Ly-1); ForceMomentsZ=zeros(Lx-1,Ly-1);
Ixx=zeros(Lx-1,Ly-1); Iyy=zeros(Lx-1,Ly-1); Izz=zeros(Lx-1,Ly-1); Ixy=zeros(Lx-1,Ly-1); Ixz=zeros(Lx-1,Ly-1); Iyz=zeros(Lx-1,Ly-1);
% Recalculate geometry to account for bending -----------------------------
[x,y,z]=RibBendingTaper(x,y,z,RibAngle,E,RibWidth,RibDepth,RibWallThickness,EndRibWidth... line 2
    ,EndRibDepth,EndRibWallThickness,SupportStrutLoc,RibLength,NRibs,nopNose,nopFrustum,nopRibLength,nopRibWidth);
%--------------------------------------------------------------------------
for i=1:Lx
    for j=1:Ly
        if i<Lx && j<Ly
            vector1=[x(i,j),y(i,j),z(i,j)];             % Define first corner of trapezoid
            vector2=[x(i+1,j),y(i+1,j),z(i+1,j)];       % Define second corner of trapezoid
            vector3=[x(i,j+1),y(i,j+1),z(i,j+1)];       % Define third corner of trapezoid
            vector4=[x(i+1,j+1),y(i+1,j+1),z(i+1,j+1)]; % Define fourth corner of trapezoid
            side1=vector2-vector1;                      % Define side vector of trapezoid
            length1=norm(side1);                        % Calculate side length of trapezoid
            side2=vector3-vector1;                      % Define top vector of trapezoid
            length2=norm(side2);                        % Calculate top length of trapezoid
            side3=vector4-vector2;                      % Define base vector of trapezoid
            length3=norm(side3);                        % Calculate base length of trapezoid
            height=sqrt(length1^2-0.25*(length3-length2)^2);    % Calculate isosceles trapezoid height
            Areas(i,j)=0.5*(length2+length3)*height;            % Calculate trapezoid area
            Volumes(i,j)=Areas(i,j)*HSthickness;                % Calculate volumes
            Masses(i,j)=Volumes(i,j)*HSdensity;                 % Calculate masses
            cross1=vector2-vector3;                     % Define diagonal vector 1 across trapezoid
            cross2=vector4-vector1;                     % Define diagonal vector 2 across trapezoid
            normals=cross(cross2,cross1);               % Calculate surface normal via cross product of diagonals
            normals=normals/norm(normals);              % Normalise surface normal
            normalx(i,j)=normals(1);                    % Place x co-ords of surface normals into matrix
            normaly(i,j)=normals(2);                    % Place y co-ords of surface normals into matrix
            normalz(i,j)=normals(3);                    % Place z co-ords of surface normals into matrix
            t1locx(i,j)=(x(i,j)+x(i+1,j)+x(i+1,j+1))/3; % Triangle 1 centroid x co-ord (for calculating trapezoid centroid)
            t1locy(i,j)=(y(i,j)+y(i+1,j)+y(i+1,j+1))/3; % Triangle 1 centroid y co-ord
            t1locz(i,j)=(z(i,j)+z(i+1,j)+z(i+1,j+1))/3; % Triangle 1 centroid z co-ord
            t2locx(i,j)=(x(i,j)+x(i,j+1)+x(i+1,j+1))/3; % Triangle 2 centroid x co-ord
            t2locy(i,j)=(y(i,j)+y(i,j+1)+y(i+1,j+1))/3; % Triangle 2 centroid y co-ord
            t2locz(i,j)=(z(i,j)+z(i,j+1)+z(i+1,j+1))/3; % Triangle 2 centroid z co-ord
            t1area=0.5*length3*height;                  % Area of triangle 1
            t2area=0.5*length2*height;                  % Area of triangle 2
            PressuresXloc(i,j)=(t1locx(i,j)*t1area+t2locx(i,j)*t2area)/Areas(i,j);   % Trapezoid centroid x co-ord
            PressuresYloc(i,j)=(t1locy(i,j)*t1area+t2locy(i,j)*t2area)/Areas(i,j);   % Trapezoid centroid y co-ord
            PressuresZloc(i,j)=(t1locz(i,j)*t1area+t2locz(i,j)*t2area)/Areas(i,j);   % Trapezoid centroid z co-ord
            %---------------this four lines take most of the time----------
            ThetasX(i,j)=atan2(norm(cross(normals,xvector)),dot(normals,xvector));    % Angle between normal and x-axis (Matlab recommended formula)
            ThetasY(i,j)=atan2(norm(cross(normals,yvector)),dot(normals,yvector));    % Angle between normal and y-axis
            ThetasZ(i,j)=atan2(norm(cross(normals,zvector)),dot(normals,zvector));    % Angle between normal and z-axis
            ThetasV(i,j)=atan2(norm(cross(normals,Vinfwind)),dot(normals,Vinfwind));  % Angle between normal and velocity vector
            %--------------------------------------------------------------
            HSAngle(i,j)=ThetasV(i,j)-pi/2;                 % Angle of heatshield to velocity vector for modified Newtonian calculation
            Cps(i,j)=Cpmax*(sin(HSAngle(i,j)))^2;           % Pressure coefficient using modified Newtonian method
            MaxPressures(i,j)=Cps(i,j)*qinf+pinf;              % Pressure calculation
            NormalForces(i,j)=MaxPressures(i,j)*Areas(i,j);          % Force normal to surface
            ForcesX(i,j)=NormalForces(i,j)*cos(ThetasX(i,j));        % Force compoment in x direction
            ForcesY(i,j)=NormalForces(i,j)*cos(ThetasY(i,j));        % Force compoment in y direction
            ForcesZ(i,j)=NormalForces(i,j)*cos(ThetasZ(i,j));        % Force component in z direction
            ForcesV(i,j)=NormalForces(i,j)*cos(ThetasV(i,j));        % Force compoment in velocity direction
            MassMomentsX(i,j)=Masses(i,j)*PressuresXloc(i,j);
            MassMomentsY(i,j)=Masses(i,j)*PressuresYloc(i,j);
            MassMomentsZ(i,j)=Masses(i,j)*PressuresZloc(i,j);
         end
    end
end
XForces=ForcesX;                      % Matrix of X forces applied to each mesh element
YForces=ForcesY;                      % Matrix of Y forces applied to each mesh element
ZForces=ForcesZ;                      % Matrix of Z forces applied to each mesh element
TotalXForces=sum(ForcesX(:));         %-TotalXForces ---> % Axial force (N)
TotalYForces=sum(ForcesY(:));         % Side force (N)
TotalZForces=sum(ForcesZ(:));         %-TotalZForces ---> % Normal force (N)
massHS=sum(Masses(:));                % Heatshield mass (kg)
massPL=mass-massHS;                                         % Payload mass (kg)
COMx=(sum(MassMomentsX(:))+massPL*COMPL(1))/mass; % X-coord of centre of mass
COMy=(sum(MassMomentsY(:))+massPL*COMPL(2))/mass; % Y-coord of centre of mass
COMz=(sum(MassMomentsZ(:))+massPL*COMPL(3))/mass; % Z-coord of centre of mass

persistent MaxXForces
if (isempty(MaxXForces)==1)
    MaxXForces = 0;   
elseif TotalXForces <= MaxXForces
    MaxXForces = TotalXForces;
    %save('MaxXForces.mat', 'MaxXForces') %takes 83s
end
    %save('MaxPressureData.mat', 'MaxPressures', 'PressuresXloc', 'PressuresYloc', 'PressuresZloc') %takes 280s
    %save('ForceData.mat','NormalForces') 
for i=1:Lx
    for j=1:Ly
        if i<Lx && j<Ly
            Ixx(i,j)=((PressuresYloc(i,j)-COMy)^2+(PressuresZloc(i,j)-COMz)^2)*Masses(i,j);
            Iyy(i,j)=((PressuresXloc(i,j)-COMx)^2+(PressuresZloc(i,j)-COMz)^2)*Masses(i,j);
            Izz(i,j)=((PressuresXloc(i,j)-COMx)^2+(PressuresYloc(i,j)-COMy)^2)*Masses(i,j);
            Ixy(i,j)=((PressuresXloc(i,j)-COMx)*(PressuresYloc(i,j)-COMy))*Masses(i,j);
            Ixz(i,j)=((PressuresXloc(i,j)-COMx)*(PressuresZloc(i,j)-COMz))*Masses(i,j);
            Iyz(i,j)=((PressuresYloc(i,j)-COMy)*(PressuresZloc(i,j)-COMz))*Masses(i,j);
            ForceMomentsX(i,j)=ForcesZ(i,j)*(PressuresYloc(i,j)-COMy)-ForcesY(i,j)*(PressuresZloc(i,j)-COMz);     % Moments of force about the Centre of Mass
            ForceMomentsY(i,j)=-ForcesZ(i,j)*(PressuresXloc(i,j)-COMx)+ForcesX(i,j)*(PressuresZloc(i,j)-COMz);    % ~Pitching moment
            ForceMomentsZ(i,j)=ForcesY(i,j)*(PressuresXloc(i,j)-COMx)-ForcesX(i,j)*(PressuresYloc(i,j)-COMy);
        end
    end
end
TotalIxx=sum(Ixx(:))+0.1535+massPL*RPayload^2/2+massPL*(COMPL(1)-COMx)^2;                          % Calculating components of moment of inertia matrix
TotalIyy=sum(Iyy(:))-0.4664+massPL/12*(3*RPayload^2+HeightPayload^2)+massPL*(COMPL(2)-COMy)^2;
TotalIzz=sum(Izz(:))+0.6236+massPL/12*(3*RPayload^2+HeightPayload^2)+massPL*(COMPL(3)-COMz)^2;
TotalIxy=sum(Ixy(:))-0.3424;
TotalIxz=sum(Ixz(:))+0.1486;
TotalIyz=sum(Iyz(:))+1.4477;
TotalFMX=sum(ForceMomentsX(:));  % Summing moments of force
TotalFMY=sum(ForceMomentsY(:));
TotalFMZ=sum(ForceMomentsZ(:));
TotalFM=[TotalFMX;TotalFMY;TotalFMZ];
IMatrix=[TotalIxx TotalIxy TotalIxz;TotalIxy TotalIyy TotalIyz;TotalIxz TotalIyz TotalIzz]; % Create moment of inertia matrix
AngAccels=IMatrix\(TotalFM-cross(omega,IMatrix*omega));              % Calculate omegadots, including cLeft matrix division \ instead of inverse
Ca=(-TotalXForces)/(qinf*RefArea);                                   % Axial force coefficient
Cn=(-TotalZForces)/(qinf*RefArea);                                   % Normal force coefficient
Cy=TotalYForces/(qinf*RefArea);                                      % Side force coefficient
Cl=Cn*cos(alpha)-Ca*sin(alpha);                                      % IPPW short course
Cd=Cn*sin(alpha)+Ca*cos(alpha);
Croll=TotalFMX/(qinf*RefArea*Diameter);                              % Rolling moment coefficient Cl
Cpitch=TotalFMY/(qinf*RefArea*Diameter);                             % Pitching moment coefficient Cm
Cyaw=TotalFMZ/(qinf*RefArea*Diameter);                               % Yawing moment coefficient Cn
LDRatio=Cl/Cd;                                                       % Lift-to-drag ratio
PsFs=[PsFs;TotalXForces,TotalYForces,TotalZForces,Cd,Cl,Ca,Cn,Cy,Croll,Cpitch,Cyaw,TotalFMX,TotalFMY,TotalFMZ,COMx,COMy,COMz];
return