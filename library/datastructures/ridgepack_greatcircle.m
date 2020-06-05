function [dist,angl,phi,tracklat,tracklon,tracklen]=...
          ridgepack_greatcircle(lat1,lon1,lat2,lon2,earthradius)

% ridgepack_greatcircle - Great circle distance and angle on sphere
%
% function [distance,angle]=...
%           ridgepack_greatcircle(lat1,lon1,lat2,lon2,earthradius)
%
% INPUT:
%
% lat1,lon2   - start point of track. This can also be entered
%               as the start of a search latitude and longitude.
%               Must be supplied in degrees north.
% lat2,lon2   - end point of track. This can be entered as a vector
%               of latitudes and longitudes to provide a vector
%               of distances. Units of degrees east.
% earthradius - Radius of the Earth. This is optional, and if
%               left out the default Ridgepack value is used.
%
% OUTPUT:
%
% dist     - distance along the great circle route in the same
%            units of earthradius between the lat-lon pairs. This 
%            has the same length has lat2 and lon2.
% angl     - arc length in degrees between the lat-lon pair.
% phi      - azimuth from (lat1,lon1) to (lat2,lon2) in degrees.
% tracklat - latitudes at 1 nautical mile intervals along track.
% tracklon - longitudes at 1 nautical mile intervals along track.
% tracklen - length along track starting from the lat1-lon1 pair
%            in the units of earthradius between the lat-lon pairs.
%            This has the same length as tracklat and tracklon
%
% Ridgepack Version 2.0
% Andrew Roberts, Los Alamos National Laboratory, 2020 
% afroberts@lanl.gov

global debug;
%debug=true;
if debug; disp(['Entering ',mfilename,'...']); end

if nargin<4
 error('not enough inputs')
elseif length(lat1)>1 | length(lon1)>1
 error('lengths of lat1 and/or lon1 exceeds 1')
elseif lat1<-90 | lat1>90
 error('Lat1 out of range')
elseif any(lat2<-90) | any(lat2>90)
 error('Lat2 out of range')
elseif lon1<-180 | lon1>180
 error('Lon1 out of range')
elseif any(lon2<-180) | any(lon2>360)
 error('Lon2 out of range')
end

if nargin<5 | isempty(earthradius)
 h=ridgepack_astroconstants;
 earthradius=h.r.const;
 if debug
  disp(['Earth Radius: ',num2str(earthradius),' ',h.r.units])
 end
end

% twist x,y,z coordinates to be centered around lat1 and lon1
[x,y,z,phi,theta]=ridgepack_satfwd(lat2,lon2,lat1,lon1,90,1,false);

% find distance along arc, and angle of arc in degrees
r=sqrt(x.^2+y.^2);
dist=earthradius.*atan2(r,z);
angl=atan2d(r,z);

% generate track points at 1 nautical mile segementations
if nargout>3
 segments=floor(60.*angl);
 if segments<=1
  tracklat=[lat1 lat2];
  tracklon=[lon1 lon2];
 else
  theta=([1:segments]/60)*(pi/180);
  z=earthradius.*cos(theta);
  x=earthradius.*sin(theta).*cos(phi);
  y=earthradius.*sin(theta).*sin(phi);
  r=sqrt(x.^2+y.^2);
  [lat,lon]=ridgepack_satinv(phi,theta,lat1,lon1);
  tracklat=[lat1 lat lat2];
  tracklon=[lon1 lon lon2];
  tracklen=[0 earthradius.*atan2(r,z) dist];
 end
end

if debug; disp(['Leaving ',mfilename,'...']); end


