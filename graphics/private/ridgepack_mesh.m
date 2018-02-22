function [latgrat,longrat,z]=ridgepack_mesh(lat,lon,z)

% ridgepack_mesh - Create mesh for mapping data from lat-lon grid
%
% function [latgrat,longrat,z]=ridgepack_mesh(lat,lon,z)
%
% This creates a grid suitable for mapping from data:
%
% Input:
% lat - latitude 1D matrix
% lon - longitude 1D matrix
% z   - array to be mapped
%
% Ouput:
% latgrat - mesh of latitude values
% longrat - mesh of longitude values
% z       - array with added values for wrapped data
%
% Andrew Roberts, Naval Postgraduate School, March 2018  (afrobert@nps.edu)
%

global debug;
if debug; disp(['Entering ',mfilename,'...']); end

if ((360.0-lon(length(lon))) == (lon(2)-lon(1)))
  if size(z,2)==length(lon)
	  z(:,length(lon)+1)=z(:,1);
  else
	  z(length(lon)+1,:)=z(1,:);
  end
  lon(length(lon)+1)=360.;
  disp('Longitudes have been wrapped');
else
  disp('Longitudes not wrapped');
end

[latgrat,longrat] = meshgrat(lat,lon);

if debug; disp(['...Leaving ',mfilename]); end

