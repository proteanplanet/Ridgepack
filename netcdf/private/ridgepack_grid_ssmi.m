function nc=ridgepack_grid_ssmi(nc,ns,n,res,xsize,ysize)

% ridgepack_grid_ssmi - Provides lat and lon for the old ssmi grid
%
% function nc=ridgepack_grid_ssmi(nc)
%
% This function provides the latitudes and longitudes
% of the old ssmi polar stereographic grid.
%
% INPUT:
%
% nc - netcdf structure to have new lats/longs and x and y data added
%      This may be left empty [] if generating a new netcdf structure.
% ns - 'north' or 'south' for north of south grids
% n  - Number or grid points to append to the boundaries of the domain
%      for the purpose of filtering data in preparation for model 
%      input fields.  This may be omitted or left empty if 
%      no grid points are to be appended to the grid.
% res- Resolution in km (optional - the default is 50km)
% xsize - Number of grid points on x direction 
% ysize - Number of grid points on y direction 
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
if nargin>1
 if isempty(n); n=0; elseif n<0; disp('NOTE: n is less than zero'); end
else
 n=0;
end

% fill in blank nc structure
if isempty(nc); nc.attributes.title='Polar Stereographic Grid'; end

if nargin<4
 deltx=50; % resolution in km
elseif isnumeric(res)
 deltx=res; % resolution in km
else
 error('res should be a number')
end

if nargin<5
 xdist=340;
elseif ~isnumeric(xsize)
 error('xsize should be a number')
else 
 xdist=xsize;
end
 
if nargin<5
 ydist=340;
elseif ~isnumeric(ysize)
 error('ysize should be a number')
else 
 ydist=ysize;
end


% PLAY WITH MAXDIST IF THERE ARE CONTOURING PROBLEMS USING ridgepack_contourM

maxdistx=xdist + 2*n;
if maxdistx<=0; 
 disp('negative n is too large'); 
elseif round(maxdistx)~=maxdistx
 error([num2str(xdist),' must be divisible by the resolution'])
end

maxdisty=ydist + 2*n;
if maxdisty<=0; 
 disp('negative n is too large');
elseif round(maxdisty)~=maxdisty
 error([num2str(ydist),' must be divisible by the resolution'])
end

disp(['Resolution is ',num2str(deltx),'km']);

if strcmp(ns,'north')
 disp('Generating ssmi grid for northern hemisphere.');
 SGN=1;
else
 disp('Generating ssmi grid for southern hemisphere.');
 SGN=-1;
end

if (nargin>5 &xsize==304 & ysize==448 & strcmp(ns,'north'))
 lonrotation=-45; % SSM/I northern 25km
 xdist=3850;
 ydist=5350;
elseif (nargin>5 &xsize==608 & ysize==896 & strcmp(ns,'north'))
 lonrotation=-45; % SSM/I northern 12.5km
 xdist=3850;
 ydist=5350;
elseif (nargin>5 &xsize==316 & ysize==332 & strcmp(ns,'south'))
 lonrotation=0; % SSM/I southern 25km
 xdist=3950;
 ydist=3950;
elseif (nargin>5 &xsize==632 & ysize==664 & strcmp(ns,'south'))
 lonrotation=0; % SSM/I southern 12.5km
 xdist=3950;
 ydist=3950;
elseif (nargin<6)
 error('Not enough arguments')
else
 lonrotation=0;
 xdist=deltx*maxdistx/2; % x distance of bottom left from South Pole
 ydist=deltx*maxdisty/2; % y distance of bottom left from South Pole
end

nc.longitude.data=-99999*ones([maxdistx maxdisty]);
nc.latitude.data=-99999*ones([maxdistx maxdisty]);

xabscissa=zeros([1 maxdistx]);
yordinate=zeros([1 maxdisty]);

for i=1:maxdistx; 
for j=1:maxdisty;

 xabscissa(i)=(i-1)*deltx-(xdist-deltx/2.);
 %yordinate(j)=(maxdisty-j)*deltx-(ydist-deltx/2.);
 yordinate(j)=(j-1)*deltx-(ydist-deltx/2.);

 [nc.latitude.data(i,j),nc.longitude.data(i,j),SLAT]=ridgepack_xytogeodetic(xabscissa(i),yordinate(j),SGN);
 
end
end

nc.latitude.data=wrapTo180(nc.latitude.data);
nc.longitude.data=wrapTo180(nc.longitude.data+lonrotation);

nc.x.data=xabscissa;
nc.x.units='km';
nc.x.long_name='x coordinate of projection';
nc.x.standard_name='projection_x_coodinate';
nc.x.dimension={'x'};
nc.x.type='NC_FLOAT';

nc.y.data=yordinate;
nc.y.units='km';
nc.y.long_name='y coordinate of projection';
nc.y.standard_name='projection_y_coodinate';
nc.y.dimension={'y'};
nc.y.type='NC_FLOAT';

nc.latitude.long_name='latitude coordinate';
nc.latitude.standard_name='latitude';
nc.latitude.dimension={'x','y'};
nc.latitude.units='degrees_north';
nc.latitude.type='NC_FLOAT';

nc.longitude.long_name='longitude coordinate';
nc.longitude.standard_name='longitude';
nc.longitude.dimension={'x','y'};
nc.longitude.units='degrees_east';
nc.longitude.type='NC_FLOAT';

nc.stereographic_mesh.grid_mapping_name='polar_stereographic';
nc.stereographic_mesh.straight_vertical_longitude_from_pole=90-SGN*90;
nc.stereographic_mesh.latitude_of_projection_origin=SGN*90;
nc.stereographic_mesh.false_easting=xdist;
nc.stereographic_mesh.false_northing=ydist;
nc.stereographic_mesh.standard_parallel=SLAT;
nc.stereographic_mesh.type='NC_CHAR';
nc.stereographic_mesh.dimension={};

nc.turn.data=SGN*(90+nc.longitude.data+lonrotation);
nc.turn.units='degrees';
nc.turn.long_name='vector turning angle to lat-long (u,v) grid from polar stereographic';
nc.turn.dimension={'x','y'};

% run check through the structure
[nc,out]=ridgepack_struct(nc);

if debug; disp(['...Leaving ',mfilename]); end

