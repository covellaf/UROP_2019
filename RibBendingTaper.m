function [x,y,z]=RibBendingTaper(x,y,z,RibAngle,E,RibWidth,RibDepth,RibWallThickness,EndRibWidth,EndRibDepth... line 2
    ,EndRibWallThickness,SupportStrutLoc,RibLength,NRibs,nopNose,nopFrustum,nopRibLength,nopRibWidth)
global count XForces YForces ZForces  MaxCurvatures MaxSlopes MaxDeltas MaxForces MaxMoments MaxMaxPrincipalStresses
% Calculate beam deflection of cantilever element of rib
indexRibs=zeros(1,NRibs); % Set up vector defining which 'js' are ribs
for j=1:NRibs+1
    indexRibs(j)=1+(j-1)*(nopRibLength-1);
end
% Calculate Bending Force mesh element strips per rib (mid-rib to mid-rib)
BendingForceMatrix=zeros(nopNose-1+nopFrustum,NRibs+1);
BendingForceXMatrix=zeros(nopNose-1+nopFrustum,NRibs+1);
BendingForceYMatrix=zeros(nopNose-1+nopFrustum,NRibs+1);
BendingForceZMatrix=zeros(nopNose-1+nopFrustum,NRibs+1);
x2=x; y2=y; z2=z;
[Li,Lj]=size(XForces);
if mod(nopRibWidth, 2) == 0         % If nopRibWidth is even (i.e. mesh elements are odd)
    for k=2:NRibs
        for j=indexRibs(k):indexRibs(k)
            for i=nopNose-1+nopFrustum:nopNose+nopFrustum+nopRibLength-3
                BendingForceXMatrix(i-nopNose-nopFrustum+2,1)=sum(XForces(i,(1:1+nopRibWidth/2-2)))+sum(XForces(i,(Lj-nopRibWidth/2+2):Lj))+XForces(i,Lj-nopRibWidth/2+1)/2+XForces(i,1+nopRibWidth/2-1)/2;
                BendingForceYMatrix(i-nopNose-nopFrustum+2,1)=sum(YForces(i,(1:1+nopRibWidth/2-2)))+sum(YForces(i,(Lj-nopRibWidth/2+2):Lj))+YForces(i,Lj-nopRibWidth/2+1)/2+YForces(i,1+nopRibWidth/2-1)/2;
                BendingForceZMatrix(i-nopNose-nopFrustum+2,1)=sum(ZForces(i,(1:1+nopRibWidth/2-2)))+sum(ZForces(i,(Lj-nopRibWidth/2+2):Lj))+ZForces(i,Lj-nopRibWidth/2+1)/2+ZForces(i,1+nopRibWidth/2-1)/2;
                BendingForceXMatrix(i-nopNose-nopFrustum+2,k)=sum(XForces(i,j+1-nopRibWidth/2:j-2+nopRibWidth/2))+XForces(i,j-nopRibWidth/2)/2+XForces(i,j-1+nopRibWidth/2)/2;
                BendingForceYMatrix(i-nopNose-nopFrustum+2,k)=sum(YForces(i,j+1-nopRibWidth/2:j-2+nopRibWidth/2))+YForces(i,j-nopRibWidth/2)/2+YForces(i,j-1+nopRibWidth/2)/2;
                BendingForceZMatrix(i-nopNose-nopFrustum+2,k)=sum(ZForces(i,j+1-nopRibWidth/2:j-2+nopRibWidth/2))+ZForces(i,j-nopRibWidth/2)/2+ZForces(i,j-1+nopRibWidth/2)/2;
            end
        end
    end
else                                % If nopRibWidth is odd (i.e. mesh elements are even)
    for k=2:NRibs
        for j=indexRibs(k):indexRibs(k)
            for i=nopNose-1+nopFrustum:nopNose+nopFrustum+nopRibLength-3
                BendingForceXMatrix(i-nopNose-nopFrustum+2,1)=sum(XForces(i,(1:1+(nopRibWidth-1)/2-1)))+sum(XForces(i,(Lj-(nopRibWidth-1)/2+1):Lj));
                BendingForceYMatrix(i-nopNose-nopFrustum+2,1)=sum(YForces(i,(1:1+(nopRibWidth-1)/2-1)))+sum(YForces(i,(Lj-(nopRibWidth-1)/2+1):Lj));
                BendingForceZMatrix(i-nopNose-nopFrustum+2,1)=sum(ZForces(i,(1:1+(nopRibWidth-1)/2-1)))+sum(ZForces(i,(Lj-(nopRibWidth-1)/2+1):Lj));
                BendingForceXMatrix(i-nopNose-nopFrustum+2,k)=sum(XForces(i,j-(nopRibWidth-1)/2:j+(nopRibWidth-1)/2-1));
                BendingForceYMatrix(i-nopNose-nopFrustum+2,k)=sum(YForces(i,j-(nopRibWidth-1)/2:j+(nopRibWidth-1)/2-1));
                BendingForceZMatrix(i-nopNose-nopFrustum+2,k)=sum(ZForces(i,j-(nopRibWidth-1)/2:j+(nopRibWidth-1)/2-1));
            end
        end
    end
end
% Sum mesh element strips into uniform force along ribs
RSSForceMatrix=sqrt(BendingForceXMatrix.^2+BendingForceYMatrix.^2+BendingForceZMatrix.^2);
FirstRow=zeros(1,size(RSSForceMatrix,2));
BendingForceXMatrix=[FirstRow;BendingForceXMatrix];
BendingForceYMatrix=[FirstRow;BendingForceYMatrix];
BendingForceZMatrix=[FirstRow;BendingForceZMatrix];
for j=1:NRibs
    Deltai=0;       % Initial deflection set to 0
    Slopei=0;       % Initial slope set to 0
    Moments(1)=count;
    MaxPrincipalStresses(1)=count;
    Curvatures(1)=count;
    Slopes(1)=count;
    Deltas(1)=count;
    Moments(2)=j;
    MaxPrincipalStresses(2)=j;
    Curvatures(2)=j;
    Slopes(2)=j;
    Deltas(2)=j;
    % Define angle matrices, bending force vectors and normal force components
    for i=nopNose-1+nopFrustum:nopNose+nopFrustum+nopRibLength-2
        Origin=[x(nopNose-1+nopFrustum,indexRibs(j)), y(nopNose-1+nopFrustum,indexRibs(j)),z(nopNose-1+nopFrustum,indexRibs(j))];   % Define origin point for this rib
        Current=[x(i,indexRibs(j)), y(i,indexRibs(j)), z(i,indexRibs(j))];
        Previous=[x(i-1,indexRibs(j)),  y(i-1,indexRibs(j)), z(i-1,indexRibs(j))];
        Between=Current-Previous;
        XAngle(i,j)=asin(Between(1)/(sqrt(Between(1)^2+Between(2)^2+Between(3)^2)));
        YAngle(i,j)=asin(Between(3)/(sqrt(Between(2)^2+Between(3)^2)));              % Calculate angles in Y-Z plane (see 07/11/17 notebook)
        ZAngle(i,j)=asin(Between(2)/(sqrt(Between(2)^2+Between(3)^2)));
        xdistCurrent(i,j)=sqrt(sum((Current - Origin) .^ 2));
        BendingForceVector=[BendingForceXMatrix(i-nopNose-nopFrustum+2,j),BendingForceYMatrix(i-nopNose-nopFrustum+2,j),BendingForceZMatrix(i-nopNose-nopFrustum+2,j)];
        NormalForce(i,j)=dot(BendingForceVector,cross(Between,cross([-1,0,0],Between))/norm(cross(Between,cross([-1,0,0],Between))));
    end
    BendingForceMatrix=NormalForce([nopNose-1+nopFrustum:end],:);
    Rssl=sum(BendingForceMatrix((1:end),j).*xdistCurrent(((nopNose-1+nopFrustum):end),j))/SupportStrutLoc;
    Rend=sum(BendingForceMatrix((1:end),j))-Rssl;
    for i=nopNose-1+nopFrustum:nopNose+nopFrustum+nopRibLength-2
        Origin=[x(nopNose-1+nopFrustum,indexRibs(j)), y(nopNose-1+nopFrustum,indexRibs(j)),z(nopNose-1+nopFrustum,indexRibs(j))];   % Define origin point for this rib
        Current=[x(i,indexRibs(j)), y(i,indexRibs(j)), z(i,indexRibs(j))];
        % Define co-ordinates of current point being evaluated
        xdistCurrent(i,j)=sqrt(sum((Current - Origin) .^ 2));                                                                             % Calculate distance between Current point and Origin point
        t=RibWallThickness-((RibWallThickness-EndRibWallThickness)*xdistCurrent(i,j)/RibLength);
        b=RibWidth-((RibWidth-EndRibWidth)*xdistCurrent(i,j)/RibLength);
        h=RibDepth-((RibDepth-EndRibDepth)*xdistCurrent(i,j)/RibLength);
        IIbeam=t*(h-2*t)^3/12+b*(h^3-(h-2*t)^3)/12;
        IRib=IIbeam;
        % Calculating initial slope required to return support strut location deflection to zero.
        if xdistCurrent(i,j)<SupportStrutLoc
            Moment=Rend*xdistCurrent(i,j)-sum(BendingForceMatrix((1:i-(nopNose-1+nopFrustum)+1),j).*(xdistCurrent(i,j)-xdistCurrent(((nopNose-1+nopFrustum):i),j)));
            MaxPrincipalStress=Moment*RibDepth/2/IRib;
            Curvature=-Moment/(E*IRib);             % Curvature, d2deltadx2
            Slopef=Slopei+Curvature*(xdistCurrent(i,j)-xdistCurrent(i-1,j));          % Slope [rad], ddeltadx
            Deltaf=Deltai+(Slopef+Slopei)/2*(xdistCurrent(i,j)-xdistCurrent(i-1,j));  % Deflection [m], delta
            Slopei=Slopef;
            Deltai=Deltaf;
            xdistFinalBeforeSS=xdistCurrent(i,j);
            k=i;
        else
            Moment=Rend*SupportStrutLoc-sum(BendingForceMatrix((1:k-(nopNose-1+nopFrustum)+1),j).*(SupportStrutLoc-xdistCurrent(((nopNose-1+nopFrustum):k),j)));
            MaxPrincipalStress=Moment*RibDepth/2/IRib;
            Curvature=-Moment/(E*IRib);             % Curvature, d2deltadx2
            Slopef=Slopei+Curvature*(SupportStrutLoc-xdistFinalBeforeSS);          % Slope [rad], ddeltadx
            Deltaf=Deltai+(Slopef+Slopei)/2*(SupportStrutLoc-xdistFinalBeforeSS);  % Deflection [m], delta
            NewSlopei=-Deltaf/SupportStrutLoc;
            break
        end
    end
    Slopei=NewSlopei;
    Deltai=0;
    % Restarting distance stepping with new initial slope to get correct deformations.
    for i=nopNose-1+nopFrustum:nopNose+nopFrustum+nopRibLength-2
        Origin=[x(nopNose-1+nopFrustum,indexRibs(j)), y(nopNose-1+nopFrustum,indexRibs(j)),z(nopNose-1+nopFrustum,indexRibs(j))];   % Define origin point for this rib
        Current=[x(i,indexRibs(j)), y(i,indexRibs(j)), z(i,indexRibs(j))];
        % Define co-ordinates of current point being evaluated
        xdistCurrent(i,j)=sqrt(sum((Current - Origin) .^ 2));                                                                             % Calculate distance between Current point and Origin point
        t=RibWallThickness-((RibWallThickness-EndRibWallThickness)*xdistCurrent(i,j)/RibLength);
        b=RibWidth-((RibWidth-EndRibWidth)*xdistCurrent(i,j)/RibLength);
        h=RibDepth-((RibDepth-EndRibDepth)*xdistCurrent(i,j)/RibLength);
        IIbeam=t*(h-2*t)^3/12+b*(h^3-(h-2*t)^3)/12;
        IRib=IIbeam;
        if xdistCurrent(i,j)<SupportStrutLoc
            Moment=Rend*xdistCurrent(i,j)-sum(BendingForceMatrix((1:i-(nopNose-1+nopFrustum)+1),j).*(xdistCurrent(i,j)-xdistCurrent(((nopNose-1+nopFrustum):i),j)));
        else
            Moment=Rend*xdistCurrent(i,j)+Rssl*(xdistCurrent(i,j)-SupportStrutLoc)-sum(BendingForceMatrix((1:i-(nopNose-1+nopFrustum)+1),j).*(xdistCurrent(i,j)-xdistCurrent(((nopNose-1+nopFrustum):i),j)));
        end
        MaxPrincipalStress=Moment*RibDepth/2/IRib;
        Curvature=-Moment/(E*IRib);             % Curvature, d2deltadx2
        Slopef=Slopei+Curvature*(xdistCurrent(i,j)-xdistCurrent(i-1,j));          % Slope [rad], ddeltadx
        Deltaf=Deltai+(Slopef+Slopei)/2*(xdistCurrent(i,j)-xdistCurrent(i-1,j));  % Deflection [m], delta
        x2(i,indexRibs(j))=x2(i,indexRibs(j))-Deltaf*sind(RibAngle);
        y2(i,indexRibs(j))=y2(i,indexRibs(j))-Deltaf*cosd(RibAngle)*sin(ZAngle(i,j));
        z2(i,indexRibs(j))=z2(i,indexRibs(j))-Deltaf*cosd(RibAngle)*sin(YAngle(i,j));
        Slopei=Slopef;
        Deltai=Deltaf;
        Moments(i)=Moment;
        MaxPrincipalStresses(i)=MaxPrincipalStress;
        Curvatures(i)=Curvature;
        Slopes(i)=Slopef;
        Deltas(i)=Deltaf;
    end
    % Extract maximum parameters
    if abs(Moments(end-2))>abs(MaxMoments(end-2))
        MaxForces=BendingForceMatrix;
        MaxMoments=Moments;
    end
    if abs(MaxPrincipalStresses(end-2))>abs(MaxMaxPrincipalStresses(end-2))
        MaxMaxPrincipalStresses=MaxPrincipalStresses;
    end
    if abs(Curvatures(end-2))>abs(MaxCurvatures(end-2))
        MaxCurvatures=Curvatures;
    end
    if abs(Slopes(end))>abs(MaxSlopes(end))
        MaxSlopes=Slopes;
    end
    if abs(Deltas(end))>abs(MaxDeltas(end))
        MaxDeltas=Deltas;
    end
end
% Match final rib to first rib
x2(:,end)=x2(:,1); y2(:,end)=y2(:,1); z2(:,end)=z2(:,1);
% Redraw geometry between ribs
t=linspace(0,1,nopRibWidth);
for k=1:NRibs
    index=indexRibs(k);
    index2=indexRibs(k+1);
    for i=1:nopRibLength
        for j=1:nopRibWidth
            indexJ=index+j-1;
            indexI=nopFrustum+nopNose+i-2;         %-2 because -2 in Rho for no redondance
            x2(indexI,indexJ)=x2(indexI,index)+t(j)*(x2(indexI,index2)-x2(indexI,index));
            y2(indexI,indexJ)=y2(indexI,index)+t(j)*(y2(indexI,index2)-y2(indexI,index));
            z2(indexI,indexJ)=z2(indexI,index)+t(j)*(z2(indexI,index2)-z2(indexI,index));
        end
    end
end
x=x2; y=y2; z=z2;
count=count+1;
return