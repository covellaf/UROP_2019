%% Defining Global variables
global x y z TrajResults PsFs MaxCurvatures MaxSlopes MaxDeltas
global XForces YForces ZForces MaxForces MaxMoments MaxMaxPrincipalStresses 
%% Importing parameters from excel file
p = readtable('GeneralData.csv.xlsx', 'ReadVariableNames',false);  
rMars = p.Var2(1);
omegaMars = p.Var2(2);
muMars = p.Var2(3);
kMars = p.Var2(4);
J2Mars = p.Var2(5);
CMars = p.Var2(6);
aMars = p.Var2(7);
bMars = p.Var2(8);
Vinf = p.Var2(9);
FPAideg = p.Var2(10);
azideg = p.Var2(11);
latdeg = p.Var2(12);
longdeg = p.Var2(13);
h = p.Var2(14);
Vorb = p.Var2(15);
Torb = p.Var2(16);
alphadeg = p.Var2(17);
betadeg = p.Var2(18);
bankdeg = p.Var2(19);
omegaXdeg = p.Var2(20);
omegaYdeg = p.Var2(21);
omegaZdeg = p.Var2(22);
Qload = p.Var2(23);
Tw = p.Var2(24);
FPA = p.Var2(25);
azi = p.Var2(26);
lat = p.Var2(27);
long = p.Var2(28);
alpha = p.Var2(29);
beta = p.Var2(30);
bank = p.Var2(31);
omegaX = p.Var2(32);
omegaY = p.Var2(33);
omegaZ = p.Var2(34);
alpha_total = p.Var2(35);
phi_a = p.Var2(36);
roll = p.Var2(37);
pitch = p.Var2(38);
yaw = p.Var2(39);           
e0 = p.Var2(40);
e1 = p.Var2(41);
e2 = p.Var2(42);
e3 = p.Var2(43);
Diameter = p.Var2(44);
RNose = p.Var2(45);
RShoulder = p.Var2(46);
SphereConeAngle = p.Var2(47);
mass = p.Var2(48);
RNoseCone = p.Var2(49);
RRib = p.Var2(50);
NRibs = p.Var2(51);
RibAngle = p.Var2(52);
RBody = p.Var2(53);
RefArea = p.Var2(54);
RPayload = p.Var2(55);
HeightPayload = p.Var2(56);
HSthickness = p.Var2(57);
RibWidth = p.Var2(58);
RibDepth = p.Var2(59);
RibWallThickness = p.Var2(60);
RootArea = p.Var2(61);
EndRibWidth = p.Var2(62);
EndRibDepth = p.Var2(63);
EndRibWallThickness = p.Var2(64);
EndArea = p.Var2(65);
RSupportStrutLoc = p.Var2(66);
RCantRib = p.Var2(67);
SupportStrutLoc = p.Var2(68);
CantRibLength = p.Var2(69);
RibLength = p.Var2(70);
E = p.Var2(71);
nu = p.Var2(72);
RibDensity = p.Var2(73);
HSdensity = p.Var2(74);
TPSspecificmass = p.Var2(75);
RootRibMass = p.Var2(76);
EndRibMass = p.Var2(77);
RibMass = p.Var2(78);
SphereConeAngleRad = p.Var2(79);                  
NoseHeight = p.Var2(80);                   
NoseLength = p.Var2(81);
ShoulderHeight = p.Var2(82);
ShoulderLength = p.Var2(83);
FrustumLength = p.Var2(84);
FrustumHeight = p.Var2(85);
ConeHeight = p.Var2(86);
ConeVertex = p.Var2(87);
RibAngleRad = p.Var2(88);
r = p.Var2(89);
u = p.Var2(90);
v = p.Var2(91);
w = p.Var2(92);
Vinfgc=[u; v; w];  
nop = p.Var2(93);
nopNose = p.Var2(94);
nopFrustum = p.Var2(95);
nopRibLength = p.Var2(96);
nopRibWidth = p.Var2(97);
nopShoulder = p.Var2(98); 
Vinfwind=[Vinf*cos(beta)*cos(alpha); Vinf*sin(beta); Vinf*cos(beta)*sin(alpha)];
COMPL=[-(0.27*Diameter); 0; 0];   
%% call Geometry and Initialisation
[x,y,z]=Geometry2(RNoseCone,RNose,RRib,NRibs,SphereConeAngleRad,RibAngleRad... line2
    ,nopNose,nopFrustum,nopRibLength,nopRibWidth,ConeVertex,NoseHeight,NoseLength);
[Lx,Ly] = size(x); %Lx = length(x)/rows and Ly = length(y)/columns
XForces=zeros(Lx-1,Ly-1);
YForces=zeros(Lx-1,Ly-1);
ZForces=zeros(Lx-1,Ly-1);
Headings={'Altitude' 'Velocity' 'MachNo' 'FlightPathAngle' 'Azimuth' 'Latitude' 'Longitude' 'DynamicPressure' 'Frustum Pressure'};
TrajResults=[[],[],[],[],[],[],[],[],[],[],[],[],[]];
PsFs=[[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[]];
MaxForces=zeros(nopNose-1+nopFrustum,NRibs+1);
MaxMoments=zeros(nop*2,1);
MaxMaxPrincipalStresses=zeros(nop*2,1);
MaxCurvatures=zeros(nop*2,1);
MaxSlopes=zeros(nop*2,1);
MaxDeltas=zeros(nop*2,1);
%% Image of heat shield and video
figure(1)
C=x; surf(x,y,z,C);
% xlabel('x'); ylabel('y');
xlabel('$x$','interpreter', 'latex', 'fontsize', 20)
ylabel('$y$','interpreter', 'latex', 'fontsize', 20)
zlabel('$z$','interpreter', 'latex', 'fontsize', 20);
axis([-10 10 -10 10 -10 10]);
colorbar
caxis([-CantRibLength 0]);
view(90,75);    % View looking slightly up on angle
global count
count=1; framecount=1;
F=getframe(gcf); imshow(F.cdata)
%% Numerical Runge-Kutta integration
tRange=[0,350];                                         % Define timespan range
var=[r;lat;long;u;v;w;e0;e1;e2;e3;omegaX;omegaY;omegaZ]; % Initial conditions for variables
    opts=odeset('Events',@EndEvents,'InitialStep',1e-8,'RelTol',1e-7); % Set break condition and initial step size
    [t,vardot]=ode45(@Trajectory,tRange,var,opts);                % Solve equations of motion
%% Fourier analysis of frequencies and generate second Image
[C,ia,ic]=unique(TrajResults(:,1),'rows');
TrajResults2=TrajResults(ia,:);
timef=TrajResults2(:,1);
anglef=TrajResults2(:,14);
timef2=(100:0.1:127)';
anglef2=interp1(timef,anglef,timef2);
anglefft=fft(anglef2);
mag=abs(anglefft);
Ts=mean(diff(timef));
Ts2=0.1;
Fs=1/Ts2;
freq=(0:(length(anglefft)-1))*Fs/length(anglefft);
figure(2)
plot(freq,mag); grid on
xlabel('$Frequency [Hz]$','interpreter', 'latex', 'fontsize', 20)
ylabel('$FFT Magnitude$','interpreter', 'latex', 'fontsize', 20)