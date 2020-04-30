function [x,y,z,phi,theta]=ridgepack_satfwd(lats,lons,centerlat,centerlon,horizon,altitude,removepoints)

% theta - polar angle
% phi - azimuthal angle

global debug;
%debug=true;
if debug; disp(['Entering ',mfilename,'...']); end

if nargin<4
 error('missing inputs')
elseif length(centerlat)>1 | length(centerlon)>1
 error('centerlat and centerlon are too long')
elseif length(lats)~=length(lons)
 error('lats and lons have unequal lengths')
end

if nargin<5 | isempty(horizon)
 horizon=90;
end

if nargin<6 | isempty(altitude)
 altitude=1;
end

if nargin<7 | ~islogical(removepoints) 
 removepoints=true;
end

% assign lats and lons
c=deg2rad(lats);
d=deg2rad(lons);

% Rodrigues' Rotation Formula
R = altitude*ones(size(c)); 
ax = R.*sin((pi/2)-c).*cos(d);
ay = R.*sin((pi/2)-c).*sin(d);
az = R.*cos((pi/2)-c);

% rotate to correct longitude
latcenter=deg2rad(90);
loncenter=deg2rad(0);
beta=deg2rad(-centerlon-90);
r=altitude;
kx=r.*sin((pi/2)-latcenter).*cos(loncenter);
ky=r.*sin((pi/2)-latcenter).*sin(loncenter);
kz=r.*cos((pi/2)-latcenter);
kdota=(kx.*ax)+(ky.*ay)+(kz.*az);
bx=cos(beta).*ax + sin(beta)*(ky*az-kz*ay) + kdota*(1-cos(beta))*kx;
by=cos(beta).*ay - sin(beta)*(kx*az-kz*ax) + kdota*(1-cos(beta))*ky;
bz=cos(beta).*az + sin(beta)*(kx*ay-ky*ax) + kdota*(1-cos(beta))*kz;

% Rotate to correct latitude
latcenter=deg2rad(0);
loncenter=deg2rad(0);
beta=deg2rad(centerlat-90);

r=altitude;
kx=r.*sin((pi/2)-latcenter).*cos(loncenter);
ky=r.*sin((pi/2)-latcenter).*sin(loncenter);
kz=r.*cos((pi/2)-latcenter);
kdota=(kx.*bx)+(ky.*by)+(kz.*bz);
cx=cos(beta).*bx + sin(beta)*(ky.*bz-kz.*by) + kdota.*(1-cos(beta))*kx;
cy=cos(beta).*by - sin(beta)*(kx.*bz-kz.*bx) + kdota.*(1-cos(beta))*ky;
cz=cos(beta).*bz + sin(beta)*(kx.*by-ky.*bx) + kdota.*(1-cos(beta))*kz;

% calculate local spherical coordinates 
phi=atan2(cy,cx);
theta=atan2(sqrt(cx.^2+cy.^2),cz);

% limit the satellite horizon
if removepoints
 idx=find(theta>deg2rad(horizon));
 phi(idx)=NaN;
 theta(idx)=NaN;
end

% calculate final x, y and z
x = r.*sin(theta).*cos(phi);
y = r.*sin(theta).*sin(phi);
z = r.*cos(theta);

if debug; disp(['Leaving ',mfilename,'...']); end

