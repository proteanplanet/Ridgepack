function [lat,lon,SLAT]=ncxytogeodetic(X,Y,SGN)

% NCXYTOGEODETIC - Convert X and Y coordinates on stereographic grid to lat-long
%
% function [lat,lon,SLAT]=ncxytogeodetic(X,Y,SGN)
% 
% This function converts x and y coordinates on a polar stereographic grid
% to latitudes and longitudes, and is essentially a JPL algorithm rewritten
% to work in Matlab.
%
% Input:
% X   - polar stereographic X coordinate (km)
% Y   - polar stereographic Y coordinate (km)
% SGN - Hemisphere specifier: 1=Northern Hemisphere, -1=Southern Hemisphere
%
% Where (X,Y)=(0,0) is at the pole.
%
% Output:
% lat   - Latitude
% lon   - Longitude
% SLAT  - Latitude of true distance
%
% *-------------------------------------------------------------------------*
% * This subroutine converts from Polar Stereographic (X,Y) coordinates     *
% * to geodetic latitude and longitude for the polar regions. The equations *
% * are from Snyder, J. P., 1982,  Map Projections Used by the U.S.         *
% * Geological Survey, Geological Survey Bulletin 1532, U.S. Government     *
% * Printing Office.  See JPL Technical Memorandum 3349-85-101 for further  *
% * details.                                                                *
% * CODE CONVERTED FROM A JPL SCRIPT BY C. S. Morris and  V. J. Troisi      *
% * With added area calculations					    *
% *-------------------------------------------------------------------------*
%
% See also: ncgeodetictoxy
%
% Converted by Andrew Roberts 2012
% Department of Oceanography, Naval Postgraduate School
%
% $Id$

global debug;
if debug; disp(['Entering ',mfilename,'...']); end

% set defaults
E2=0.006693883;
E=sqrt(E2); % eccentricity of Hughes ellipsoid
RE=6378.273; % earth radius in km
SLAT=70; % latitude of true distance
PI=pi; % PI
CDR=180/pi; %  Conversion constant from degrees to radians

SL = SLAT*PI/180;

RHO=sqrt((X.^2)+(Y.^2));

CM=cos(SL)/sqrt(1-E2*(sin(SL).^2));
T=tan((PI/4)-(SL/2))./(((1.0-E*sin(SL))./(1.0+E*sin(SL))).^(E/2));

if abs(SLAT-90.)<1.E-5
 T=(RHO.*sqrt(((1+E).^(1+E)).*((1-E).^(1-E)))/2)./RE; 
else
 T=(RHO*T)./(RE*CM);
end

CHI=(PI/2)-(2*atan(T));
ALAT=CHI+((E2/2.0)+(5.0*E2.^2.0/24.0)+(E2.^3.0/12.0))*sin(2*CHI)+...
            ((7.0*E2.^2.0/48.0)+(29.0*E2.^3/240.0))*sin(4.0*CHI)+...
            (7.0*E2.^3.0/120.0)*sin(6.0*CHI);
ALAT=SGN*ALAT;
ALONG=atan2(SGN*X,-SGN*Y);
ALONG=SGN*ALONG;

ALAT(RHO<0.1)=pi*SGN/2;
ALONG(RHO<0.1)=0.0;

%-------------------------------------------------------------------------*

lon=wrapTo180(ALONG*180./pi);
lat=wrapTo180(ALAT*180./pi);

if debug; disp(['...Leaving ',mfilename]); end
