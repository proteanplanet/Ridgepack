function [z,nc]=ridgepack_polarplot(nc,Z)

% ridgepack_polarplot - Regrids data to a polar stereographic grid for polar plots
%
% function [z,nc]=ridgepack_polarplot(nc,Z)
%
% This function regrids data to a polar stereographic grid for viewing 
% if a polar stereographic projection is used, and the data is on 
% 1D lat and 1D long described grid (lat-long grid).
%
% INPUT:
%
% nc  - netcdf structure
% Z   - name of variable to be extracted
%
%
% OUTPUT:
%
% z   - data to be interpolated
% nc  - netcdf structure re-gridded
%
% Ridgepack Version 1.0
% Andrew Roberts, Naval Postgraduate School, March 2018 (afrobert@nps.edu)

global debug;
if debug; disp(['Entering ',mfilename,'...']); end

% get map axis information
h=getm(gca);

% check to see if this is stereographic, and if we are on a lat long grid
if (size(nc.latitude.data,1)==1 | size(nc.latitude.data,2)==1) && strcmp(h.mapprojection,'stereo');
        if h.origin(1)>70;
                nc=ridgepack_regrid(nc,Z,'',2);
                z=nc.(Z).data;
                disp('Using northern SSM/I grid for stereographic projection'); 
        elseif h.origin(1)<-70;
                nc=ridgepack_regrid(nc,Z,'',3);
                z=nc.(Z).data;
                disp('Using southern SSM/I grid for stereographic projection'); 
        else
                if debug; disp('Not using SSM/I grid for stereographic projection'); end
        end
else
	if strcmp(h.mapprojection,'stereo')
		if debug; disp('Not regridding to SSM/I grid for stereographic projection'); end
	end
	z=nc.(Z).data;
end

if debug; disp(['...Leaving ',mfilename]); end

