function ridgepack_sunlightm(time)

% ridgepack_sunlightm - Generates a camera light on a 3D globe map plot
%
% function ridgepack_sunlightm(time)
%
% Adds sun light to a  3D plot of earth based on the time
% provided through the global time variable.
%
% Input:
% time - Serial time in UTC.  This may be omitted for the time to
% be set to 0300UTC on the Boreal Summer Solstice. Otherwise the 
% solar angle is roughly calculated for the time provided.
% maph - axis handle for globe map.
%
% Andrew Roberts, Naval Postgraduate School, March 2018  (afrobert@nps.edu)
%

global debug;
if debug; disp(['Entering ',mfilename,'...']); end

if ~ismap(gca)
 error('Must be applied to a current map handle')
end

h=gcm;

if ~all(strcmp(h.mapprojection,'globe')) ;
 error('ridgepack_sunlightm set only for globe plots');
end

lighting none;

if ~exist('time') || isempty(time)
 sunangle=26;
 sunlongitude=-90;
else 
 disp(['Sunlight for: ',datestr(time)]);
 [Y, M, D, H, MN, S] = datevec(time);
 if Y<1800 ; error('Year is too early to correctly calulate sun angle'); end
 doy=time-datenum(Y, 1, 1, 0, 0, 0);
 sunangle=-23.45*cosd((doy+10)*360/365);
 sunlongitude=180-360*(time-floor(time));
end

lightm(sunangle,sunlongitude, 0);
material dull;
lighting gouraud;

if debug; disp(['...Leaving ',mfilename]); end



