function [distance,angle,tracklat,tracklon,length]=...
          ridgepack_greatcircle(lat1,lon1,lat2,lon2,earthradius)

% ridgepack_greatcircle - Great circle distance and angle on sphere
%
% function [distance,angle]=...
%           ridgepack_greatcircle(lat1,lon1,lat2,lon2,earthradius)
%
% INPUT:
%
% lat1,lon2   - start point of track
% lat2,lon2   - end point of track
% earthradius - Radius of the Earth. This is optional, and if
%               left out the default Ridgepack value is used.
%
% OUTPUT:
%
% distance - distance along the great circle route in the same
%            units of earthradius between the lat-lon pairs.
% angle    - arc length in degrees between the lat-lon pair.
% tracklat - latitudes at 1 nautical mile intervals along track.
% tracklon - longitudes at 1 nautical mile intervals along track.
% length   - length along track starting from the lat1-lon1 pair
%            in the units of earthradius between the lat-lon pairs.
%            This has the same length as tracklat and tracklon

if nargin<5 | isempty(earthradius)
 h=ridgepack_astroconstants;
 earthradius=h.r.const;
 disp(['Earth Radius Units: ',h.r.units])
elseif nargin<4
 error('not enough inputs')
elseif lat1<-90 | lat1>90
 error('Lat1 out of range')
elseif lat2<-90 | lat2>90
 error('Lat2 out of range')
elseif lon1<-180 | lon1>180
 error('Lon1 out of range')
elseif lon2<-180 | lon2>180
 error('Lon2 out of range')
end

% twist x,y,z coordinates to be centered around lat1 and lon1
[x,y,z,phi,theta]=ridgepack_satfwd(lat2,lon2,lat1,lon1,90,1,false);

% find distance along arc, and angle of arc in degrees
r=sqrt(x.^2+y.^2);
distance=earthradius*atan2(r,z);
angle=atan2d(r,z);

% generate track points at 1 nautical mile segementations
if nargout>2
 segments=floor(60*angle);
 if segments<=1
  tracklat=[lat1 lat2];
  tracklon=[lon1 lon2];
 else
  theta=([1:segments]/60)*(pi/180);
  z=earthradius*cos(theta);
  x=earthradius*sin(theta)*cos(phi);
  y=earthradius*sin(theta)*sin(phi);
  r=sqrt(x.^2+y.^2);
  [lat,lon]=ridgepack_satinv(phi,theta,lat1,lon1,90,1);
  tracklat=[lat1 lat lat2];
  tracklon=[lon1 lon lon2];
  length=[0 earthradius*atan2(r,z) distance];
 end
end

