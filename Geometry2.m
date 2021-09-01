function [x,y,z]=Geometry2(RNoseCone,RNose,RRib,NRibs,SphereConeAngleRad,RibAngleRad... line2
    ,nopNose,nopFrustum,nopRibLength,nopRibWidth,ConeVertex,NoseHeight,NoseLength)

RShoulder=0; nopShoulder=0;
ShoulderHeight=RShoulder*(1-cos(SphereConeAngleRad));    
ShoulderLength=RShoulder*sin(SphereConeAngleRad);
FrustumLength=RNoseCone-NoseLength-ShoulderHeight;
FrustumHeight=FrustumLength/tan(SphereConeAngleRad);
r1=linspace(0,NoseLength,nopNose);                    % Define radial vector for nose
r2=linspace(NoseLength,RNoseCone,nopFrustum);         % Define radial vector for frustum
r3=linspace(RNoseCone,RNoseCone+RRib,nopRibLength);   % Define radial vector for cloth
if ShoulderHeight==0
    r4=[];
else
    r4=linspace(RNoseCone+RRib,RNoseCone+RRib+ShoulderHeight,nopShoulder); % Define radial vector for shoulder (if included)
end
Rho=[r1,r2(2:end),r3(2:end),r4(2:end)];                         % Combine to create single vector for full heatshield radius
nopRho = length(Rho);                                           % Number of mesh points for radial length
nopTheta = NRibs*(nopRibWidth-1);                               % Number of mesh points for revolution
Theta=linspace(0,2*pi,nopTheta+1);                              % Define theta angle vector for nosecone heatshield
[theta,rho]=meshgrid(Theta,Rho);                                % Create mesh grid in r and theta co-ordinates
z2=zeros(nopRho,nopTheta+1);                                    % Set up empty matrix
for i=1:nopRho
    for j=1:nopTheta+1
        if rho(i,j)<=NoseLength
            z2(i,j)=-sqrt(abs(RNose^2-rho(i,j).^2))+RNose;
        elseif rho(i,j)>NoseLength && rho(i,j)<=(FrustumLength+NoseLength)
            z2(i,j)=rho(i,j)/tan(SphereConeAngleRad)+ConeVertex;
        elseif rho(i,j)>FrustumLength+NoseLength && rho(i,j)<=(FrustumLength+NoseLength+RRib)
            z2(i,j)=rho(i,j)/tan(RibAngleRad) -( Rho(nopFrustum+nopNose)/tan(RibAngleRad)- z2(nopFrustum+nopNose-1,j))+ (r3(2)-r3(1))/tan(RibAngleRad);
        elseif rho(i,j)<=RNoseCone
            z2(i,j)=-sqrt(abs(RShoulder^2-(rho(i,j)-(RNoseCone-RShoulder)).^2))+(ShoulderLength+FrustumHeight+NoseHeight);
        end
    end
end
indexRibs=zeros(1,NRibs);                   % Set up vector defining which j points are ribs
for j=1:NRibs+1                             % +1 because it is the same as the first one 1
    indexRibs(j)=1+(j-1)*(nopRibWidth-1);
end
x=-z2;
y=rho.*cos(theta);
z=rho.*sin(theta);
t=linspace(0,1,nopRibWidth);
for k=1:NRibs
    index=indexRibs(k);
    index2=indexRibs(k+1);   
    for i=1:nopRibLength 
        for j=1:nopRibWidth
            indexJ=index+j-1;
            indexI=nopFrustum+nopNose+i-2;         %-2 because -2 in Rho for no redondance
            x(indexI,indexJ)=x(indexI,index)+t(j)*(x(indexI,index2)-x(indexI,index));
            y(indexI,indexJ)=y(indexI,index)+t(j)*(y(indexI,index2)-y(indexI,index));
            z(indexI,indexJ)=z(indexI,index)+t(j)*(z(indexI,index2)-z(indexI,index));
        end
    end
end
end