function [lat,lon,h,k]=ridgepack_xytoease(s,r,s0,r0,C,SGN)

% ridgepack_xytoease - Convert indices on Azimuthal Polar EASE grid to lat-long
%
% function [lat,lon,h,k]=ridgepack_xytoease(r,s,SGN)
% 
% This function converts x and y coordinates on the EASE grid to latitudes 
% and longitudes, where EASE stands for Equal-Area Scalable Earth. For 
% more information, see http://nsidc.org/data/ease/ease_grid.html.
% This function is for Azimuthal Polar EASE grids.
%
% INPUT:
%
% s   - vector of EASE row coordinates
% r   - vector of EASE column coordinates 
% s0  -	Map origin row
% r0  -	Map origin column
% C   -	Nominal cell (pixel) size in km
% SGN - Hemisphere specifier: 1=Northern Hemisphere, -1=Southern Hemisphere
%
%
% OUTPUT:
%
% lat   - Latitude
% lon   - Longitude
% h     - scale along meridian
% k     - scale along parallel
%
% Note:  
% For several datasets available at the National Snow and Ice Data Center (NSIDC) 
% with documented resolution of 25km on the EASE grid, the 'real' resolution, as 
% input into the NSIDC published EASE grid formulas is 25.067525 km. For more 
% information, see http://nsidc.org/data/pathfinders/grid.html. 
%
% Converted by Andrew Roberts 2013
%

global debug;
if debug; disp(['Entering ',mfilename,'...']); end

% set defaults
R=6371.228; % Radius of the Earth (km)

for i=1:length(s)
for j=1:length(r); 

 lambda(i,j)=atan2((r(i)-r0),(s(j)-s0));
 if abs(sin(lambda(i,j)))<eps 
   phi(i,j)=(pi/2)-2*asin(C*(s(j)-s0)/(2*R*cos(lambda(i,j))));
 else
   phi(i,j)=(pi/2)-2*asin(C*(r(i)-r0)/(2*R*sin(lambda(i,j))));
 end
 h(i,j)=cos(pi/4-phi(i,j)/2);
 k(i,j)=sec(pi/4-phi(i,j)/2);

 lambda(i,j)=SGN*lambda(i,j)+pi/2;

 lon(i,j)=wrapTo180(lambda(i,j)*180/pi);
 lat(i,j)=wrapTo180(SGN*phi(i,j)*180/pi);

end; 
end;

if debug; disp(['...Leaving ',mfilename]); end
