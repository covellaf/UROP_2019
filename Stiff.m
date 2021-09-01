function EI=Stiff(j,xdistCurrent,E,RibWallThickness,EndRibWallThickness,RibWidth,EndRibWidth,RibDepth,EndRibDepth,RibLength,option,const)
%% The purpose of the function is to give the stiffness of the element of beam
% j                       -   integer  - Index of the rib
% EI                      -   float    - Stiffness of the element considered
% RibWallThickness        -   vector(1xNRibs)- Thickness of the ribs
% EndRibWallThickness     -   vector(1xNRibs)- Thickness of the ribs
% RibWidth                -   vector(1xNRibs)- Width of the ribs
% EndRibWidth             -   vector(1xNRibs)- Width of the ribs
% RibDepth                -   vector(1xNRibs)- Depth of the ribs
% EndRibDepth             -   vector(1xNRibs)- Depth of the ribs
% RibLength               -   vector(1xNRibs)- Length of the ribs
% xdistCurrent            -   float     - Position of the lement considered
% option                  -   interger  - 1: Stiff Constant equal to Const; 2: Stiff computed
%% Using the given constant stiffness EI=const
if option==1
    EI=const;
%% Computing EI(j,z) thanks to the shape of the beams and the young modulus
elseif option==2
    t=RibWallThickness-((RibWallThickness-EndRibWallThickness).*xdistCurrent./RibLength);
    b=RibWidth-((RibWidth-EndRibWidth).*xdistCurrent./RibLength);
    h=RibDepth-((RibDepth-EndRibDepth).*xdistCurrent./RibLength);
    IRib=t.*(h-t.*2).^3./12+b.*(h.^3-(h-t.*2).^3)./12;
    %IRib(j)=0.0017;
    %ISolid=b*h^3/12;
    %IHollow=b*h^3/12-((b-2*t)*(h-2*t)^3/12);
    %ISolidTube=pi/4*(h/2)^4;
    %ITube=pi/4*((h/2)^4-(h/2-t)^4;
    %ITbeam=b*t^3/12
    %ILbeam=t*(5*h^2-5*h*t+t^2)*(h^2-h*t+t^2)/(12*(2*h-t));
    EI=E(j)*IRib(j);
end
%% Function in progress(UROP Student), to compute the stiffness given EI at the origin/end plus the order of the expression
I1
I2
'linear'
if variation=='linear'
I(position)=(I2-I1)*xdistcurrent/RibLength+I1
EI=E*I
elseif variation=='quadratic'
    
end
