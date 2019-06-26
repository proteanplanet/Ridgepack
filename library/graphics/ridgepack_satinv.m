function [lats,lons]=ridgepack_satinv(phi,theta,centerlat,centerlon,horizon,altitude)

% Rodrigues' Rotation Formula
R = altitude*ones(size(phi)); 
ax = R.*sin(theta).*cos(phi);
ay = R.*sin(theta).*sin(phi);
az = R.*cos(theta);

% Rotate back to correct latitude
latcenter=deg2rad(0);
loncenter=deg2rad(0);
beta=deg2rad(centerlat-90);
r=altitude;
kx=r.*sin((pi/2)-latcenter).*cos(loncenter);
ky=r.*sin((pi/2)-latcenter).*sin(loncenter);
kz=r.*cos((pi/2)-latcenter);
kdota=(kx.*ax)+(ky.*ay)+(kz.*az);
cx=cos(beta).*ax+sin(beta)*(ky.*az-kz.*ay)+kdota.*(1-cos(beta))*kx;
cy=cos(beta).*ay-sin(beta)*(kx.*az-kz.*ax)+kdota.*(1-cos(beta))*ky;
cz=cos(beta).*az+sin(beta)*(kx.*ay-ky.*ax)+kdota.*(1-cos(beta))*kz;

% rotate back to correct longitude
latcenter=deg2rad(90);
loncenter=deg2rad(0);
beta=deg2rad(-centerlon+90);
r=altitude;
kx=r.*sin((pi/2)-latcenter).*cos(loncenter);
ky=r.*sin((pi/2)-latcenter).*sin(loncenter);
kz=r.*cos((pi/2)-latcenter);
kdota=(kx.*cx)+(ky.*cy)+(kz.*cz);
bx=cos(beta).*cx+sin(beta)*(ky*cz-kz*cy)+kdota*(1-cos(beta))*kx;
by=cos(beta).*cy-sin(beta)*(kx*cz-kz*cx)+kdota*(1-cos(beta))*ky;
bz=cos(beta).*cz+sin(beta)*(kx*cy-ky*cx)+kdota*(1-cos(beta))*kz;

% calculate local spherical coordinates 
phi=atan2(by,bx);
theta=atan2(sqrt(bx.^2+by.^2),bz);

% calculate lats and lons 
lats=rad2deg(theta+(pi/2));
lons=rad2deg(phi);

