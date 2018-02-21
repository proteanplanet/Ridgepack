function nc=ncgrid_rasmwrf(nc,n)

% NCGRID_RASMWRF - Provides lat and lon for the RASM WRF grid
%
% function nc=ncgrid_rasmwrf(nc,n)
%
% This function provides the latitudes and longitudes
% of the RASM WRF model.
%
% Input:
% nc - netcdf structure to have new lats/longs and x and y data added
%      This may be left empty [] if generating a new netcdf structure.
% n  - number of grid points to append along each boundary of the grid 
%      for the purpose of filtering data in preparation for model 
%      input fields.  This may be omitted or left empty if 
%      no grid points are to be appended to the grid.
%
% Output:
% nc - netcdf structure with new lats/longs and x and y added
% 
% Written by Andrew Roberts 2012
% Department of Oceanography, Naval Postgraduate School
%
% $Id$

global debug;
if debug; disp(['Entering ',mfilename,'...']); end

disp('Generating RASM WRF grid');

% check input
if nargin>1
  error('Cannot currently specify n outside of model bounds for WRF grid')
% if isempty(n); n=0; elseif n<0; error('n must be greater than 0'); end
else
 n=0;
end

disp('******CHECK TURNING ANGLE*****');

ncnew=ncclone('~/data/MODEL/RACM/WRF_CELL_AREA.nc');

nc.x=ncnew.x;
nc.y=ncnew.y;
nc.longitude=ncnew.longitude;
nc.latitude=ncnew.latitude;
nc.cell_area=ncnew.cell_area;

nc.turn.data=-wrapTo180(nc.longitude.data+114);
nc.turn.units='degrees';
nc.turn.long_name='vector turning angle to lat-long (u,v) grid from ice-tide grid';
nc.turn.type='NC_FLOAT';
nc.turn.dimension={'y','x'};

if debug; disp(['...Leaving ',mfilename]); end

