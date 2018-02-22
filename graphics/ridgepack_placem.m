function [lat,lon]=ridgepack_placem(Name)

% ridgepack_placem - Returns the lat and long of a named city
%
% function [lat,lon]=ridgepack_placem(Name)
%
% Returns the latitude and longitude of the named place
%
% Input:
% Name - city name
%
% Output:
% lat - latitude
% lon - longitude
%
% Andrew Roberts, Naval Postgraduate School, March 2018  (afrobert@nps.edu)
%

global debug;
if debug; disp(['Entering ',mfilename,'...']); end

cities = shaperead('worldcities.shp', 'UseGeoCoords', true);
numb = strmatch(Name,{cities(:).Name});

if not(isempty(numb)) ;
        lat=cities(numb).Lat;
        lon=cities(numb).Lon;
elseif strmatch(Name,'Hobart');
        lat=-42.86;
        lon=147.32;
else
        disp('No cities found for that name');
        return
end

disp(['latitude=',num2str(lat),'   longitude=',num2str(lon)]);

drawnow

if debug; disp(['...Leaving ',mfilename]); end


