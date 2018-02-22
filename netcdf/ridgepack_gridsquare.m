function [nc,xdist,ydist]=ridgepack_gridsquare(nc,model)

% ridgepack_gridsquare - Generate vertices from generalized curvilinear grid cell central grid points
%
% function [nc,xdist,ydist]=ridgepack_gridsquare(nc,model)
%
% This function calculates the vertices for the rectangular tessallation of 
% a latitude-longitude grid, thus effectively calculating the 'split grid'
% for the B-grid center or U-V points. 
%
% Input:
% nc    - netcdf structure (see ridgepack_struct for more information).
%         This structure must be in generalized lat-lon coordinates
%         as described in the help page for ridgepack_sph2gen output. An example
%         of the appropriate input is provided here:
% netcdf structure nc {
%  attributes {global}
%     Conventions: 'CF-1.0'
%     history: '2007-11-14 22:43:55 GMT by mars2netcdf-0.92'
%      
%  6 variables, dimension=[ ], data=( )
%     
%  time [01-Sep-1957 12:00:00 to 01-Aug-2002 12:00:00, timestep 31 days]
%     long_name: 'time'
%         units: 'hours since 1900-01-01 00:00:0.0'
%    dimension: {'time'}
%         data: [540x1 double]
%         type: 'NC_FLOAT'
%                                                 
%  x [1 to 144, range 143, step 1]
%    long_name: 'x-coordinate in Cartesian system'
%    dimension: {'x'}
%         data: [1x144 double]
%         type: 'NC_FLOAT'
%                    
%  y [1 to 73, range 72, step 1]
%    long_name: 'y-coordinate in Cartesian system'
%    dimension: {'y'}
%         data: [1x73 double]
%         type: 'NC_FLOAT'
%                                                               
%  latitude (-90 to 90, range 180, delta 2.5)
%    long_name: 'latitude'
%        units: 'degrees_north'
%    dimension: {'y'  'x'}
%         data: [73x144 double]
%         type: 'NC_FLOAT'
%
%  longitude (0 to 357.5, range 357.5, delta 2.5)
%    long_name: 'longitude'
%        units: 'degrees_east'
%    dimension: {'y'  'x'}
%         data: [73x144 double]
%         type: 'NC_FLOAT'
%                                                                     
%  p2t (198.5297 to 321.1956, range 122.6659, median 283.0118)
%    long_name: '2 metre temperature'
%        units: 'K'
%    dimension: {'time'  'y'  'x'}
%         data: [540x73x144 double]
%    fillvalue: -32767
%         type: 'NC_FLOAT'
%  coordinates: 'latitude longitude'
%
% }
%
% model - type of model to use to calculate area on the grid integer:
%         1: Assumes the earth is a sphere
%         2: Assumes the wgs84 ellipsoid
%         3: Assumes a polar stereographic plane
%
% Output:
% nc - netcdf structure with a split grid, or corner coordinates
%      for a grid with latitude and longitude coordinates. An example
%      is provided here.
%
% netcdf structure nc {
%  attributes {global}
%     Conventions: 'CF-1.0'
%     history: '2007-11-14 22:43:55 GMT by mars2netcdf-0.92'
%	      
% 10 variables, dimension=[ ], data=( )
%     
% time [01-Sep-1957 12:00:00 to 01-Aug-2002 12:00:00, timestep 31 days]
%     long_name: 'time'
%         units: 'hours since 1900-01-01 00:00:0.0'
%     dimension: {'time'}
%          data: [540x1 double]
%          type: 'NC_FLOAT'
%		  
%  x [1 to 144, range 143, step 1]
%     long_name: 'x-coordinate in Cartesian system'
%     dimension: {'x'}
%          data: [1x144 double]
%          type: 'NC_FLOAT'
%				     
%  x_corner [0.5 to 144.5, range 144, step 1]
%     long_name: 'x-coordinate vertex in Cartesian system'
%     dimension: {'x_corner'}
%          data: [1x145 double]
%          type: 'NC_FLOAT'
%								        
%  y [1 to 73, range 72, step 1]
%     long_name: 'y-coordinate in Cartesian system'
%     dimension: {'y'}
%          data: [1x73 double]
%          type: 'NC_FLOAT'
%									   
%  y_corner [0.5 to 73.5, range 73, step 1]
%     long_name: 'y-coordinate vertex in Cartesian system'
%     dimension: {'y_corner'}
%          data: [1x74 double]
%          type: 'NC_FLOAT'
%												      
%  latitude (-90 to 90, range 180, delta 2.5)
%     long_name: 'latitude'
%         units: 'degrees_north'
%     dimension: {'y'  'x'}
%          data: [73x144 double]
%          type: 'NC_FLOAT'
%								 
%  latitude_corner (-90 to 90, range 180, median 0)
%     long_name: 'latitude of corner grid points'
%         units: 'degrees_north'
%     dimension: {'y_corner'  'x_corner'}
%          data: [74x145 double]
%          type: 'NC_FLOAT'
%
%  longitude (0 to 357.5, range 357.5, delta 2.5)
%     long_name: 'longitude'
%         units: 'degrees_east'
%     dimension: {'y'  'x'}
%          data: [73x144 double]
%          type: 'NC_FLOAT'
%
%  longitude_corner (0.33095 to 359.0811, range 358.7502, median 178.75)
%     long_name: 'longitude of corner grid points'
%         units: 'degrees_east'
%     dimension: {'y_corner'  'x_corner'}
%          data: [74x145 double]
%          type: 'NC_FLOAT'
%												  
%  p2t (198.5297 to 321.1956, range 122.6659, median 283.0118)
%     long_name: '2 metre temperature'
%         units: 'K'
%     dimension: {'time'  'y'  'x'}
%          data: [540x73x144 double]
%     fillvalue: -32767
%          type: 'NC_FLOAT'
%   coordinates: 'latitude longitude'
%
% }
%
% Andrew Roberts, Naval Postgraduate School, March 2018  (afrobert@nps.edu)
%

global debug;
if debug; disp(['Entering ',mfilename,'...']); end

% Check that nc is a structure
if not(isstruct(nc));
 error([inputname(1),' is not a structure']);
end

if nargin<2
 model=1;
elseif ~isnumeric(model)
 error('model should be a number')
end

% run a check on the data and make sure corner and middle arrays are aligned
if (isfield(nc,'x_corner') & isfield(nc,'x')) | ...
   (isfield(nc,'y_corner') & isfield(nc,'y')) | ...
   (isfield(nc,'latitude_corner') & isfield(nc,'latitude')) | ...
   (isfield(nc,'longitude_corner') & isfield(nc,'longitude'))
 if debug; disp('A split grid appears to be in this nc structure'); end
 nc=ridgepack_shuffle(nc,{'y_corner','x_corner'});
 nc=ridgepack_shuffle(nc,{'y','x'});
 return;
else
 disp('Generating a split geolocated grid from the existing grid')
end

% Run quickly through ridgepack_sph2gen to check for correct setup
nc=ridgepack_sph2gen(nc);

% make sure there is an x and y field defining rectangular grid cells
if ~isfield(nc,'x') | ~isfield(nc,'y')
 error('x and/or y are not fields')
elseif ~(min(size(nc.x.data))==1) | ~(min(size(nc.y.data))==1)
 error('x or y data are multidimensional')
elseif ~strcmp('x',nc.x.dimension) | ~strcmp('y',nc.y.dimension)
 error('The dimensions of x and/or y are incorrect')
end

% determine if the grid is global or regional
gridtype=ridgepack_gridtype(nc);

% shuffle to put x last
nc=ridgepack_shuffle(nc,{'x','y'});

% check sizes of two grids are the same and get the dimensions
sizex=size(nc.latitude.data);
sizey=size(nc.longitude.data);
if sizex~=sizey
 error(['Grid sizes differ between latitude and longitude']);
end
dimi=sizex(find(strcmp(nc.latitude.dimension,'x')));
dimj=sizex(find(strcmp(nc.latitude.dimension,'y')));
if max(size(nc.x.data))~=dimi
 disp(['Size of x is ',num2str(size(nc.x.data)),' and size of lat or lon is ',num2str(dimi)])
 error('Dimension of x coordinate does not agree with x side of lat and lon')
elseif max(size(nc.y.data))~=dimj
 disp(['Size of y is ',num2str(size(nc.y.data)),' and size of lat or lon is ',num2str(dimj)])
 error('Dimension of y coordinate does not agree with y side of lat and lon')
end

% first case of global rectangular grid in curvilinear coorinates
if strcmp('global',gridtype)

 if debug; disp('Processing global grid'); end

 % establish dimensions
 x_corner=0.5:1:dimi+0.5;
 y_corner=0.5:1:dimj+0.5;
 lon_corner=ones(dimi+1, dimj+1);
 lat_corner=ones(dimi+1, dimj+1);

 % generate large array containing boundary conditions identical to next-in value
 lon=ones(dimi+2,dimj+2);
 lat=ones(dimi+2,dimj+2);

 lon(1:dimi,1:dimj)=nc.longitude.data(:,:);
 lat(1:dimi,1:dimj)=nc.latitude.data(:,:);

 lon(3:dimi+2,1:dimj)=nc.longitude.data(:,:);
 lat(3:dimi+2,1:dimj)=nc.latitude.data(:,:);

 lon(1:dimi,3:dimj+2)=nc.longitude.data(:,:);
 lat(1:dimi,3:dimj+2)=nc.latitude.data(:,:);

 lon(3:dimi+2,3:dimj+2)=nc.longitude.data(:,:);
 lat(3:dimi+2,3:dimj+2)=nc.latitude.data(:,:);

 lon(2:dimi+1,2:dimj+1)=nc.longitude.data(:,:);
 lat(2:dimi+1,2:dimj+1)=nc.latitude.data(:,:);

 % determine along which dimension the zipline exists and the pole for global dataset
 xdist=distance(lat(:,1),lon(:,1),lat(:,end),lon(:,end),earthRadius('km'));
 ydist=distance(lat(1,:),lon(1,:),lat(end,:),lon(end,:),earthRadius('km'));

 % now generate wrapped boundary conditions along stich line and at pole
 if median(xdist)<median(ydist); 

  lat(:,1)=lat(:,dimj+1);
  lon(:,1)=lon(:,dimj+1);
  lat(:,end)=lat(:,2);
  lon(:,end)=lon(:,2);

  if model==1
   [plat1,plon1]=meanm(nc.latitude.data(1,:),nc.longitude.data(1,:),earthRadius('km'));
  elseif model==2
   [plat1,plon1]=meanm(nc.latitude.data(1,:),nc.longitude.data(1,:),wgs84Ellipsoid('km'));
  else
   error('No Earth Model Radius provided for this value of model')
  end

  if debug; disp(['Lower pole is at: ',num2str([plat1,plon1])]); end
  lat(1,2:dimj+1)=plat1;
  lon(1,2:dimj+1)=plon1;
  lat_corner(1,:)=plat1;

  if model==1
   [plat2,plon2]=meanm(nc.latitude.data(end,:),nc.longitude.data(end,:),earthRadius('km'));
  elseif model==2
   [plat2,plon2]=meanm(nc.latitude.data(end,:),nc.longitude.data(end,:),wgs84Ellipsoid('km'));
  else
   error('No Earth Model Radius provided for this value of model')
  end

  if debug; disp(['Upper pole is at: ',num2str([plat2,plon2])]); end
  lat(end,2:dimj+1)=plat2;
  lon(end,2:dimj+1)=plon2;
  lat_corner(end,:)=plat2;

 else

  lat(1,:)=lat(dimi+1,:);
  lon(1,:)=lon(dimi+1,:);
  lat(end,:)=lat(2,:);
  lon(end,:)=lon(2,:);

  if model==1
   [plat1,plon1]=meanm(nc.latitude.data(:,1),nc.longitude.data(:,1),earthRadius('km'));
  elseif model==2
   [plat1,plon1]=meanm(nc.latitude.data(:,1),nc.longitude.data(:,1),wgs84Ellipsoid('km'));
  else
   error('No Earth Model Radius provided for this value of model')
  end

  if debug; disp(['Lower pole is at: ',num2str([plat1,plon1])]); end
  lat(2:dimi+1,1)=plat1;
  lon(2:dimi+1,1)=plon1;
  lat_corner(:,1)=plat1;

  if model==1
   [plat2,plon2]=meanm(nc.latitude.data(:,end),nc.longitude.data(:,end),earthRadius('km'));
  elseif model==2
   [plat2,plon2]=meanm(nc.latitude.data(:,end),nc.longitude.data(:,end),wgs84Ellipsoid('km'));
  else
   error('No Earth Model Radius provided for this value of model')
  end

  if debug; disp(['Upper pole is at: ',num2str([plat2,plon2])]); end
  lat(2:dimi+1,end)=plat2;
  lon(2:dimi+1,end)=plon2;
  lat_corner(:,end)=plat2;

 end

 % now stack the lats and longs for each four grid corners on top
 % of each other in a column to precondition before calculating the
 % lats and longs of the corner grid points using meanm
 newlat=zeros([4 dimi+1 dimj+1]);
 newlon=zeros([4 dimi+1 dimj+1]);

 newlat(1,:,:)=lat(1:end-1,1:end-1);
 newlat(2,:,:)=lat(1:end-1,2:end);
 newlat(3,:,:)=lat(2:end,1:end-1);
 newlat(4,:,:)=lat(2:end,2:end);

 newlon(1,:,:)=lon(1:end-1,1:end-1);
 newlon(2,:,:)=lon(1:end-1,2:end);
 newlon(3,:,:)=lon(2:end,1:end-1);
 newlon(4,:,:)=lon(2:end,2:end);


 % calculate the corner lats and longs
 if model==1
  [lat_corner,lon_corner]=meanm(newlat(:,:),newlon(:,:),earthRadius('km'));
 elseif model==2
  [lat_corner,lon_corner]=meanm(newlat(:,:),newlon(:,:),wgs84Ellipsoid('km'));
 else
  error('No Earth Model Radius provided for this value of model')
 end


 % reshape lat_corner and lon_corner from 1D meanm output
 lat_corner=reshape(lat_corner,dimi+1,dimj+1);
 lon_corner=reshape(lon_corner,dimi+1,dimj+1);

% second case of limited area grid
elseif strcmp('regional',gridtype)

 if debug; disp('Processing regional grid'); end

 % establish dimensions and base arrays
 x_corner=ones([dimi+1 1]);
 y_corner=ones([dimj+1 1]);
 lon_corner=ones(dimi+1, dimj+1);
 lat_corner=ones(dimi+1, dimj+1);


 % populate x_corner and y_corner
 x_corner(1)=nc.x.data(1)-(nc.x.data(2)-nc.x.data(1))/2;
 y_corner(1)=nc.y.data(1)-(nc.y.data(2)-nc.y.data(1))/2;
 x_corner(end)=nc.x.data(end)+(nc.x.data(end)-nc.x.data(end-1))/2;
 y_corner(end)=nc.y.data(end)+(nc.y.data(end)-nc.y.data(end-1))/2;
 x_corner(2:end-1)=(nc.x.data(1:end-1)+nc.x.data(2:end))/2;
 y_corner(2:end-1)=(nc.y.data(1:end-1)+nc.y.data(2:end))/2;


 % check the projection and set the SGN
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


 % generate geographic coordinates
 if model==3

  if SGN==-99999;
   SGN=sign(mean(nc.latitude.data(:)));
  end
 
  disp('Calculating latitudes and longitudes on polar stereographic plane')

  for i=1:dimi+1
  for j=1:dimj+1
   [lat_corner(i,j),lon_corner(i,j),SLAT]=ridgepack_xytogeodetic(x_corner(i),y_corner(j),SGN);
  end
  end

  % apply grid rotation angle of original grid
  [lat,lon]=ridgepack_xytogeodetic(nc.x.data(1),nc.y.data(1),SGN);
  lon_corner=lon_corner-(lon-nc.longitude.data(1,1));

 else

  % generate large lat/lon array containing boundary conditions
  lon=ones(dimi+2, dimj+2);
  lat=ones(dimi+2, dimj+2);

  % bottom-left boundary vertex
  lon(1,1)=nc.longitude.data(1,1);
  lat(1,1)=nc.latitude.data(1,1);

  % top-left boundary vertex
  lon(dimi+2,1)=nc.longitude.data(end,1);
  lat(dimi+2,1)=nc.latitude.data(end,1);

  % bottom-right boundary vertex
  lon(1,dimj+2)=nc.longitude.data(1,end);
  lat(1,dimj+2)=nc.latitude.data(1,end);

  % top-right boundary vertex
  lon(dimi+2,dimj+2)=nc.longitude.data(end,end);
  lat(dimi+2,dimj+2)=nc.latitude.data(end,end);

  % left boundary side
  lon(2:dimi+1,1)=nc.longitude.data(:,1);
  lat(2:dimi+1,1)=nc.latitude.data(:,1);

  % right boundary side
  lon(2:dimi+1,end)=nc.longitude.data(:,end);
  lat(2:dimi+1,end)=nc.latitude.data(:,end);

  % bottom boundary side
  lon(1,2:dimj+1)=nc.longitude.data(1,:);
  lat(1,2:dimj+1)=nc.latitude.data(1,:);

  % top boundary side
  lon(end,2:dimj+1)=nc.longitude.data(end,:);
  lat(end,2:dimj+1)=nc.latitude.data(end,:);

  % now fill in array center 
  lon(2:dimi+1,2:dimj+1)=nc.longitude.data(:,:);
  lat(2:dimi+1,2:dimj+1)=nc.latitude.data(:,:);

  % now stack the lats and longs for each four grid corners on top
  % of each other in a column to precondition before calculating the
  % lats and longs of the corner grid points using meanm

  newlat=zeros([4 dimi+1 dimj+1]);
  newlon=zeros([4 dimi+1 dimj+1]);

  newlat(1,:,:)=lat(1:end-1,1:end-1);
  newlat(2,:,:)=lat(1:end-1,2:end);
  newlat(3,:,:)=lat(2:end,1:end-1);
  newlat(4,:,:)=lat(2:end,2:end);

  newlon(1,:,:)=lon(1:end-1,1:end-1);
  newlon(2,:,:)=lon(1:end-1,2:end);
  newlon(3,:,:)=lon(2:end,1:end-1);
  newlon(4,:,:)=lon(2:end,2:end);

  % calculate the corner lats and longs
  if model==1
   [lat_corner,lon_corner]=meanm(newlat(:,:),newlon(:,:),earthRadius('km'));
  elseif model==2
   [lat_corner,lon_corner]=meanm(newlat(:,:),newlon(:,:),wgs84Ellipsoid('km'));
  else
   error('wrong earth model chosen')
  end

  % reshape lat_corner and long_corner from 1D meanm output
  lat_corner=reshape(lat_corner,dimi+1,dimj+1);
  lon_corner=reshape(lon_corner,dimi+1,dimj+1);

 end

else

 error([gridtype,' not recognized for creating a split grid']);

end

% make all longitudes 0 to 360;
lat_corner=wrapTo180(lat_corner);
lon_corner=wrapTo360(lon_corner);

% add the grid to the netcdf structure
if isfield(nc,'x_corner') | isfield(nc,'y_corner') 
 error('x_corner and/or y_corner already variables in the netcdf structure')
end
nc=ridgepack_add(nc,'x_corner',x_corner','x-coordinate vertex in Cartesian system',{'x_corner'},'','NC_FLOAT');
nc=ridgepack_add(nc,'y_corner',y_corner','y-coordinate vertex in Cartesian system',{'y_corner'},'','NC_FLOAT');
if isfield(nc,'latitude_corner') | isfield(nc,'longitude_corner') 
 error('latitude_corner and/or longitude_corner already variables in the netcdf structure')
end
nc=ridgepack_add(nc,'latitude_corner',lat_corner,'latitude of corner grid points',{'x_corner','y_corner'},'degrees_north');
nc=ridgepack_add(nc,'longitude_corner',lon_corner,'longitude of corner grid points',{'x_corner','y_corner'},'degrees_east');
nc.latitude_corner.type=nc.latitude.type;
nc.longitude_corner.type=nc.longitude.type;


% run check through the structure
[nc,out]=ridgepack_struct(nc);

if debug; disp(['...Leaving ',mfilename]); end

