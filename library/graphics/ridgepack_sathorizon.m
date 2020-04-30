function [h]=ridgepack_sathorizon(centlat,centlon,horizon,...
                                  satlat,satlon,sathorizon,color)

% ridgepack_sathorizon - Draw a satellite horizon on a map
%
% function [h]=ridgepack_sathorizon(centlat,centlon,horizon,...
%
%                                  satlat,satlon,sathorizon,color)
% 
% Can be used for creating search horizon

% generate horizon in local coordinates
phi=deg2rad([0:0.001:361]);
theta=deg2rad(sathorizon*ones(size(phi)));
altitude=1;

% generate lat and longs from local satellite horizon
[lats,lons]=ridgepack_satinv(phi,theta,satlat,satlon,sathorizon,altitude);

% generate coordinates on main satellite view of globe
[x,y,z]=ridgepack_satfwd(lats,lons,centlat,centlon,90,altitude);

% plot the horizon
h=plot3(x,y,1.01*z,'Color',color)

