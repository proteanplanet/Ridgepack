function [X,Y,SLAT,K]=ridgepack_geodetictoxy(lat,lon,SGN)

% ridgepack_geodetictoxy - Convert lat-long to X and Y coordinates on a stereographic grid
%
% function [X,Y,SLAT,K]=ridgepack_geodetictoxy(lat,lon,SGN)
%
% This function converts latitudes and longitudes to x and y coordinates
% on a polar stereographic projection, and is essentially a JPL algorithm
% rewritten to work in Matlab.
% 
% Input:
% lat  - latitude (degrees North)
% lon  - longitude (degrees East)
% SGN  - Hemisphere specifier: 1=Northern Hemisphere, -1=Southern Hemisphere
%
% Output:
% X   - polar stereographic X coordinate (km)
% Y   - polar stereographic Y coordinate (km)
% SLAT - latitude of true distance (default is 70 degrees if omitted)
% K   - length scale on the stereographic projection
%
% Where (X,Y)=(0,0) is at the pole.
%
% *-------------------------------------------------------------------------*
% * This subroutine converts from geodetic latitude and longitude to Polar  *
% * Stereographic (X,Y) coordinates for the polar regions.  The equations   *
% * are from Snyder, J. P., 1982,  Map Projections Used by the U.S.         *
% * Geological Survey, Geological Survey Bulletin 1532, U.S. Government     *
% * Printing Office.  See JPL Technical Memorandum 3349-85-101 for further  *
% * details.                     			                    *
% * CODE CONVERTED FROM A JPL SCRIPT BY C. S. Morris and  V. J. Troisi      *
% *-------------------------------------------------------------------------*
%
% In addition to Morris and Troisi routine, a calculation of scale has also been 
% made using equations on pp 161 and 315 of Snyder, J.P., 1987, "Map Projections: 
% A Working Manual", USGS. Note that the equation for calculating K at the pole
% on page 161 (equation 21-35) incorrectly includes the radius of the Earth.
%
% See also: ridgepack_xytogeodetic
%
% Converted by Andrew Roberts 2012
%

global debug;
if debug; disp(['Entering ',mfilename,'...']); end

% set defaults
E2=0.006693883;
E=sqrt(E2); % eccentricity of Hughes ellipsoid
RE=6378.273; % earth radius in km
SLAT=70; % latitude of true distance
PI=pi; % PI
CDR=180/pi; %  Conversion constant from degrees to radians

% convert latitude and longitude to radians
lat=abs(lat)*PI./180;
lon=lon*PI./180;

T=tan(PI/4-lat/2)./((1-E*sin(lat))./(1+E*sin(lat))).^(E/2);
M=cos(lat)./sqrt(1-E2*(sin(lat).^2));

SL=SLAT*PI/180;
TC=tan(PI/4-SL/2)./((1-E*sin(SL))./(1+E*sin(SL))).^(E/2);
MC=cos(SL)/sqrt(1-E2*(sin(SL).^2));

if abs(90-SLAT)<1.E-5;
 RHO=2*RE*T./sqrt(((1+E)^(1+E))*((1-E)^(1-E)));
else
 RHO=RE*MC*T./TC;
end

Y=-SGN*RHO.*cos(SGN*lon);
X= SGN*RHO.*sin(SGN*lon);

X(abs(lat)>=PI/2)=0.0;
Y(abs(lat)>=PI/2)=0.0;

if RHO==0; % Calculate scale at the pole
 K=0.5*(MC./TC)*sqrt(((1+E)^(1+E))*((1-E)^(1-E)));
else % or elsewhere
 K=RHO./(RE*M);
end

if debug; disp(['...Leaving ',mfilename]); end

