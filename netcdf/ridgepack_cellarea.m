function [nc]=ridgepack_cellarea(nc,model)

% ridgepack_cellarea - Esimates the area grid cells for model or observational datasets
%
% function [nc]=ridgepack_cellarea(nc,model)
%
% This function calculates the area of each grid cell on the sphere, wgs84 
% ellipsoid or polar stereographic plane for a geolocated grid, either for 
% a global mesh, or limited area domain.  The cell area is added to the 
% netcdf structure, which also recognizes generalized coordinates from a 
% lat-long grid if need be (see ridgepack_sph2gen for more information). The output 
% is the same spherical or generalized coordinates that was input in nc 
% (see ridgepack_sph2gen for an explanation of each coordinate structure).
%
% Input:
% nc - netcdf structure (see ridgepack_struct for more information) with 
%      geolocated data.
% model - type of model to use to calculate area on the grid integer:
%      1: Assumes the earth as a sphere
%      2: Assumes the wgs84 ellipsoid
%      3: Assumes a polar stereographic plane
%
% Output:
% nc - netcdf structure with a cell_area added to the structure. This is
%      an element that looks like this in the netcdf structure:
%
% ...
% ...
% cell_area (3.7647e-07 to 0.00015163, range 0.00015125, median 0.00010711)
%       long_name: 'area on the sphere/wgs84 ellipsoid'
%       dimension: {'y'  'x'}
%            data: [73x144 double]
%            type: 'NC_DOUBLE'
%     coordinates: 'latitude longitude'
% ...
% ...
% 
% Note that the calculations are done by first calculating a split grid from the 
% vertices of four-vertex grid cells, and then assuming great circles between 
% these vertices to calculate the area within each grid cell provided in the 
% latitude and longitude arrays in the nc structure input.
%
% Andrew Roberts, Naval Postgraduate School, March 2018  (afrobert@nps.edu)
%

global debug;
if debug; disp(['Entering ',mfilename,'...']); end

disp('Calculating the area on the sphere of each grid point.');

% check inputs
if not(isstruct(nc));
 error([inputname(1),' is not a structure']);
end

if nargin<2
 model=1;
elseif ~isnumeric(model)
 error('model should be a number')
end

% check that latitude and longitude are variables
if ~isfield(nc,'longitude') | ~isfield(nc,'latitude')
 disp('No latitude & longitude coordinates to calculate cell area')
 return
elseif isfield(nc,'cell_area')
 disp('A cell_area variable already exists for this netcdf structure')
 return
end

% Convert to generalized coordinates if need be
[nc,calc]=ridgepack_sph2gen(nc);

% create a split grid from the existing grid 
ncsquare=ridgepack_gridsquare(nc,model);

% check for corner values
if ~isfield(ncsquare,'latitude_corner') | ~isfield(ncsquare,'longitude_corner')
 error('This function set up for four vertices per grid cell only');
end

% check sizes of arrays 
sizex=size(nc.latitude.data);
sizey=size(nc.longitude.data);
if sizex~=sizey
 error(['Grid sizes differ between ',x,' and ',y]);
end
dimj=sizex(1);
dimi=sizex(2);
if max(size(nc.x.data))~=dimi
  error('Dimension of x coordinate does not agree with x side of lat and lon')
elseif max(size(nc.y.data))~=dimj
  error('Dimension of y coordinate does not agree with x side of lat and lon')
end

SGN=-99999;
[variablenames,numbervariables]=ridgepack_name(nc);
for i=1:numbervariables
 name=char(variablenames{i});
 if isfield(nc.(name),'grid_mapping_name') && ...
    strcmp(nc.(name).grid_mapping_name,'polar_stereographic');
      if isfield(nc.(name),'latitude_of_projection_origin')
       SGN=sign(nc.stereographic_mesh.latitude_of_projection_origin);
      end
      if model~=3
       disp(['NOTE: Grid mapping is ',nc.(name).grid_mapping_name])
       disp('Switching to model=3 to geolocate grid edges')
       model=3;
      end
 end
end

if model==1 | model==2

 % assign the lat and long corner grid points to local arrays
 lat_corner=ncsquare.latitude_corner.data;
 lon_corner=ncsquare.longitude_corner.data;

 % calculate area on the sphere of each rectangular tile of a curvilinear grid
 nanarray=nan([dimj dimi]);
 llat=zeros([6 dimj dimi]);
 llon=zeros([6 dimj dimi]);
 llat(1,:,:)=lat_corner(1:dimj,1:dimi);
 llat(2,:,:)=lat_corner(1:dimj,2:dimi+1);
 llat(3,:,:)=lat_corner(2:dimj+1,2:dimi+1);
 llat(4,:,:)=lat_corner(2:dimj+1,1:dimi);
 llat(5,:,:)=lat_corner(1:dimj,1:dimi);
 llat(6,:,:)=nanarray(1:dimj,1:dimi);
 llon(1,:,:)=lon_corner(1:dimj,1:dimi);
 llon(2,:,:)=lon_corner(1:dimj,2:dimi+1);
 llon(3,:,:)=lon_corner(2:dimj+1,2:dimi+1);
 llon(4,:,:)=lon_corner(2:dimj+1,1:dimi);
 llon(5,:,:)=lon_corner(1:dimj,1:dimi);
 llon(6,:,:)=nanarray(1:dimj,1:dimi);

 if model==1
  disp('Using Earth as a sphere')
  cell_area=areaint(llat(:,:),llon(:,:),earthRadius('km'));
  cell_area=reshape(cell_area,[dimj dimi]);
 elseif model==2
  disp('Using WGS 84 ellipsoid')
  cell_area=areaint(llat(:,:),llon(:,:),wgs84Ellipsoid('km'));
  cell_area=reshape(cell_area,[dimj dimi]);
 end

elseif model==3

  disp('Using Polar Stereographic values to calculate areas')

  if SGN==-99999;
   SGN=sign(mean(nc.latitude.data(:)));
  end

  cell_area=zeros([dimj dimi]);

  for i=1:dimi
  for j=1:dimj
   [X,Y,SLAT,K]=ridgepack_geodetictoxy(nc.latitude.data(j,i),nc.longitude.data(j,i),SGN);
   cell_area(j,i)=abs(ncsquare.x_corner.data(i+1)-ncsquare.x_corner.data(i))*...
                  abs(ncsquare.y_corner.data(j+1)-ncsquare.y_corner.data(j))/(K^2);
  end
  end

else

  error(['No Earth Model Radius provided for this value of model=',num2str(model)])

end


% now assign the value to the nc grid
if model==1
 nc=ridgepack_add(nc,'cell_area',cell_area,'area of grid cell on sphere',nc.longitude.dimension,'km^2');
elseif model==2
 nc=ridgepack_add(nc,'cell_area',cell_area,'area of grid cell on wgs84 ellipsoid',nc.longitude.dimension,'km^2');
elseif model==3
 nc=ridgepack_add(nc,'cell_area',cell_area,'area of grid cell on stereographic plane',nc.longitude.dimension,'km^2');
else
 error('model choice does not apply')
end
nc.cell_area.standard_name='area';
nc.cell_area.coordinates='latitude longitude';

% change back to spherical coordinates if need be
if calc; nc=ridgepack_gen2sph(nc); end

[nc,out]=ridgepack_struct(nc);

disp('Finished calculating the area on the sphere of each grid point.');

if debug; disp(['...Leaving ',mfilename]); end


