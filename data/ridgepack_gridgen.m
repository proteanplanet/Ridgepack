function [nc]=ridgepack_gridgen(nc,mode,n,res)

% ridgepack_gridgen - Generate a preset grid in a netcdf structure
%
% function [nc]=ridgepack_gridgen(nc,mode,n,res)
%
% This function generates a pre-specified geographic grid
% and places the lats, longs, x and y into the netcdf
% output structure nc.
%
% INPUT:
%
% nc   - nc structure. This may be left empty [] if generating 
%        a new netcdf structure based on the grid.
%
% mode - mode=1: ice-tide grid
%        mode=2: Square polar stereographic grid, north pole, 50km resolution 
%        mode=3: Square polar stereographic grid, south pole, 50km resolution
%        mode=4: Square polar stereographic grid, north pole, 5km resolution
%        mode=5: Square polar stereographic grid, variable resolution
%        mode=6: RASM-WRF grid
%        mode=7: RASM-POP/CICE grid
%	 mode=8; SMMR & SSM/I 25km grid, north pole
%	 mode=9; SMMR & SSM/I 25km grid, south pole
%	 mode=10; AMSR-E 12.5km grid, north pole
%	 mode=11; AMSR-E 12.5km grid, south pole
%	 mode=12; 1/10th degree global spherical grid
%
% n    - number of grid points to append to the boundaries
%        of the grid.  Does not work for mode=8 to 11
%
% res  - explicitly set resolution in km for SSM/I grids
%        if res is negative, it chooses SH, else NH is 
%        the default.  Does not work for mode=8 to 11
%
% Ridgepack Version 1.0
% Andrew Roberts, Naval Postgraduate School, March 2018 (afrobert@nps.edu)
%

global debug;
if debug; disp(['Entering ',mfilename,'...']); end

if isempty(nc)
 out=true;
else
 out=false;
end

% get new lats and longs of grid
if mode==1
 nc=ridgepack_grid_icetide(nc,n);
 if out; nc.attributes.title='Hibler/Roberts Model Grid'; end
elseif mode==2
 nc=ridgepack_grid_ssmi(nc,'north',n,50,340,340);
 if out; nc.attributes.title='Polar Stereographic Northern Grid'; end
 nc=ridgepack_cellarea(nc,3);
elseif mode==3
 nc=ridgepack_grid_ssmi(nc,'south',n,50,340,340);
 if out; nc.attributes.title='Polar Stereographic Southern Grid'; end
 nc=ridgepack_cellarea(nc,3);
elseif mode==4
 nc=ridgepack_grid_ssmi(nc,'north',n,5,3400,3400);
 if out; nc.attributes.title='Polar Stereographic Northern Grid'; end
 nc=ridgepack_cellarea(nc,3);
elseif mode==5
 if res>0
  nc=ridgepack_grid_ssmi(nc,'north',n,abs(res));
  if out; nc.attributes.title='Polar Stereographic Northern Grid'; end
  nc=ridgepack_cellarea(nc,3);
 else
  nc=ridgepack_grid_ssmi(nc,'south',n,abs(res));
  if out; nc.attributes.title='Polar Stereographic Southern Grid'; end
  nc=ridgepack_cellarea(nc,3);
 end
elseif mode==6
 nc=ridgepack_grid_rasmwrf(nc);
 if out; nc.attributes.title='RASM WRF Grid'; end
elseif mode==7
 nc=ridgepack_grid_rasmpop(nc);
 if out; nc.attributes.title='RASM POP/CICE Grid'; end
elseif mode==8
 nc=ridgepack_grid_ssmi(nc,'north',0,25,304,448);
 if out; nc.attributes.title='SMMR & SSM/I Polar Stereographic Northern Grid'; end
 nc=ridgepack_cellarea(nc,3);
elseif mode==9
 nc=ridgepack_grid_ssmi(nc,'south',0,25,316,332);
 if out; nc.attributes.title='SMMR & SSM/I Polar Stereographic Southern Grid'; end
 nc=ridgepack_cellarea(nc,3);
elseif mode==10
 nc=ridgepack_grid_ssmi(nc,'north',0,12.5,608,896);
 if out; nc.attributes.title='AMSR-E Polar Stereographic Northern Grid'; end
 nc=ridgepack_cellarea(nc,3);
elseif mode==11
 nc=ridgepack_grid_ssmi(nc,'south',0,12.5,632,664);
 if out; nc.attributes.title='AMSR-E Polar Stereographic Southern Grid'; end
 nc=ridgepack_cellarea(nc,3);
elseif mode==12
 nc=ridgepack_grid_gaussian(nc,10);
 if out; nc.attributes.title='1/10th Degree Spherical Grid'; end
else
 error('No grid available for that mode number');
end

if out; ridgepack_struct(nc); end

if debug; disp(['...Leaving ',mfilename]); end

