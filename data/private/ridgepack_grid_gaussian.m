function nc=ridgepack_grid_gaussian(nc,resolution)

% ridgepack_grid_gaussian - Provides lat and lon for the old ssmi grid
%
% function nc=ridgepack_grid_gaussian(nc,resolution)
%
% This function provides the latitudes and longitudes
% of a global spherical grid.
%
% INPUT:
%
% nc - netcdf structure to have new lats/longs and x and y data added
%      This may be left empty [] if generating a new netcdf structure.
% resolution - resolution of grid in terms of number of grid points
%      per degree of the spherical grid.
%
%
% OUTPUT:
%
% nc - netcdf structure with new lats/longs and x and y added
% 
% Ridgepack Version 1.0
% Andrew Roberts, Naval Postgraduate School, March 2018 (afrobert@nps.edu)
%

global debug;
if debug; disp(['Entering ',mfilename,'...']); end

% check input
if nargin==1
 error('missing spherical grid resolution')
elseif nargin==0
 error('missing input')
elseif nargin==2 & ~isnumeric(resolution)
 error('resolution must be a numerical value')
end

% fill in blank nc structure
if isempty(nc); nc.attributes.title='Global Spherical Grid'; end

disp(['Resolution of global spherical grid is ',num2str(1/resolution),' degree']);

R=[resolution 90 -180];
Z=ones(resolution*[180 360]);

nc.x.data=[1:resolution*360];
nc.x.long_name='x-coordinate of project';
nc.x.dimension={'x'};
nc.x.type='NC_FLOAT';

nc.y.data=[1:resolution*180];
nc.y.long_name='y coordinate of projection';
nc.y.dimension={'y'};
nc.y.type='NC_FLOAT';

nc.latitude.long_name='latitude coordinate';
nc.latitude.standard_name='latitude';
nc.latitude.dimension={'y','x'};
nc.latitude.units='degrees_north';
nc.latitude.type='NC_FLOAT';

nc.longitude.long_name='longitude coordinate';
nc.longitude.standard_name='longitude';
nc.longitude.dimension={'y','x'};
nc.longitude.units='degrees_east';
nc.longitude.type='NC_FLOAT';

[nc.latitude.data nc.longitude.data]=meshgrat(Z,R);

% run check through the structure
[nc,out]=ridgepack_struct(nc);

if debug; disp(['...Leaving ',mfilename]); end

