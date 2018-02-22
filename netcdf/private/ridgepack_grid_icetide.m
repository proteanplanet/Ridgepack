function nc=ridgepack_grid_icetide(nc,n)

% ridgepack_grid_icetide - Provides lat and lon for the ice-tide grid
%
% function nc=ridgepack_grid_icetide(nc)
%
% This function provides the latitudes and longitudes
% of the model used for Bill HIbler's ice-tide work.
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
% Andrew Roberts, Naval Postgraduate School, March 2018  (afrobert@nps.edu)
%

global debug;
if debug; disp(['Entering ',mfilename,'...']); end

disp('Generating icepack grid');

% check input
if nargin>1
 if isempty(n); n=0; elseif n<0; error('n must be greater than 0'); end
else
 n=0;
end

dx=13.89; % resolution of square grid cells
nx=440+2*n; % distance in grid points along x axis
ny=292+2*n; % distance in grid points along y axis
fe=180+n; % distance in grid point from bottom of x axis
fn=160+n; % distance in grid point from top of y axis

for i=1:ny; 
for k=1:nx;

	ax=k-fe;ay=fn-i;
	r=sqrt(ax*ax+ay*ay);
	yy=57.3*atan2(ay,ax)+40.;
	if (yy<0) ; 
		yy=yy+360. ; 
        end;
	if (yy>360) ; 
		yy=yy-360. ; 
	end;
	fi(k,ny-i+1)=90.-r/8.;
	dl(k,ny-i+1)=yy;
	if (ax==0 & ay==0); 
		dl(k,ny-i+1)=0.; 
		fi(k,ny-i+1)=90. ; 
	end;
end
end

nc.x.data=dx*[0:1:nx-1]';
nc.x.units='km';
nc.x.long_name='x coordinate of projection';
nc.x.standard_name='projection_x_coodinate';
nc.x.dimension={'x'};
nc.x.type='NC_FLOAT';

nc.y.data=dx*[0:1:ny-1]';
nc.y.units='km';
nc.y.long_name='y coordinate of projection';
nc.y.standard_name='projection_y_coodinate';
nc.y.dimension={'y'};
nc.y.type='NC_FLOAT';

nc.latitude.data=fi';
nc.latitude.long_name='latitude coordinate';
nc.latitude.standard_name='latitude';
nc.latitude.units='degrees_north';
nc.latitude.type='NC_FLOAT';
nc.latitude.dimension={'y','x'};

nc.longitude.data=dl';
nc.longitude.long_name='longitude coordinate';
nc.longitude.standard_name='longitude';
nc.longitude.units='degrees_east';
nc.longitude.type='NC_FLOAT';
nc.longitude.dimension={'y','x'};

nc.icepack_mesh.grid_mapping_name='azimuthal_equidistant';
nc.icepack_mesh.longitude_of_projection_origin=-50;
nc.icepack_mesh.latitude_of_projection_origin=90;
nc.icepack_mesh.false_easting=fe*dx;
nc.icepack_mesh.false_northing=(ny-fn)*dx;
nc.icepack_mesh.type='NC_CHAR';
nc.icepack_mesh.dimension={};

nc.turn.data=(90+nc.longitude.data-50);
nc.turn.units='degrees';
nc.turn.long_name='vector turning angle to lat-long (u,v) grid from ice-tide grid';
nc.turn.type='NC_FLOAT';
nc.turn.dimension={'y','x'};

if debug; disp(['...Leaving ',mfilename]); end

