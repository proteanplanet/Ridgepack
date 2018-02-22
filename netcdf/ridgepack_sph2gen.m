function [nc,calc]=ridgepack_sph2gen(nc)

% ridgepack_sph2gen - Convert nc structure from lat-lon dimensions to general coordinates
%
% function [nc,calc]=ridgepack_sph2gen(nc)
%
% This function adds 2D latitude and longitude variables to 
% a netcdf structure, placing the existing latitude and longitude data
% into 2D arrays. The reason this is useful is that it allows routines 
% that reduce the dimensions of netcdf data and manipulate the data to treat 
% all linear geospatial data in the same manner within ncbox. 
%
% INPUT:
%
% nc - netcdf structure (see ridgepack_struct for more information) with data in 
%      spherical coordinates. An examples of the lat-lon dimension structure
%      is provided here:
%
% netcdf structure nc {
%  attributes {global}
%     Conventions: 'CF-1.0'
%         history: '2007-11-14 22:43:55 GMT by mars2netcdf-0.92'
%	      
%  4 variables, dimension=[ ], data=( )
%	       
%  latitude [-90 to 90, range 180, step -2.5]
%    long_name: 'latitude'
%        units: 'degrees_north'
%    dimension: {'latitude'}
%         data: [73x1 double]
%         type: 'NC_FLOAT'
%					  
%  longitude [0 to 357.5, range 357.5, step 2.5]
%    long_name: 'longitude'
%        units: 'degrees_east'
%    dimension: {'longitude'}
%         data: [144x1 double]
%         type: 'NC_FLOAT'
%						     
%  time [01-Sep-1957 12:00:00 to 01-Aug-2002 12:00:00, timestep 31 days]
%    long_name: 'time'
%        units: 'hours since 1900-01-01 00:00:0.0'
%    dimension: {'time'}
%         data: [540x1 double]
%         type: 'NC_FLOAT'
%				        
%  p2t (198.5297 to 321.1956, range 122.6659, median 283.0118)
%    long_name: '2 metre temperature'
%        units: 'K'
%    dimension: {'time'  'latitude'  'longitude'}
%         data: [540x73x144 double]
%    fillvalue: -32767
%         type: 'NC_FLOAT'
%  
% }
%
%
% OUTPUT:
%
% nc   - netcdf structure (see ridgepack_struct for more information) with data 
%        placed in a standardised generalised coordinates on a sphere.  
%        An example of the structure created from the above input structure 
%        is provided below.
%
% calc - logical telling if calculations were performed on the structure (calc=true)
%        or whether the structure was already in generalized coordinates (calc=false).
%
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
%
% Ridgepack Version 1.0
% Andrew Roberts, Naval Postgraduate School, March 2018 (afrobert@nps.edu)
%

global debug;
if debug; disp(['Entering ',mfilename,'...']); end

% Check that nc is a structure
if not(isstruct(nc));
 error([inputname(1),' is not a structure']);
end

% assume that calculations are to be done.
calc=true;

% Check that nc is in spherical coordinates
if ~isfield(nc,'longitude') | ~isfield(nc,'latitude')
 error('This netcdf structure is not in geographical coordinates')
elseif length(nc.longitude.dimension)==1 & ...
       strcmp(nc.longitude.dimension,char(nc.latitude.dimension)) & ...
       length(nc.latitude.dimension)==1 & ...
       strcmp(nc.latitude.dimension,char(nc.longitude.dimension)) 
  disp('Latitude and longitude describe a series')
  return
elseif length(nc.longitude.dimension)>1 & ~any(strcmp('longitude',nc.longitude.dimension)) & ...
       length(nc.latitude.dimension)>1 & ~any(strcmp('latitude',nc.latitude.dimension)) 
 if(debug);disp('This nc structure already has 2D lat and long arrays');end
 calc=false;
 if any(strcmp(nc.latitude.dimension,'x')) && any(strcmp(nc.latitude.dimension,'y'))
  return
 elseif any(strcmp(nc.latitude.dimension,'time')) | any(strcmp(nc.longitude.dimension,'time'))
  return
 else
  name{1}=char(nc.latitude.dimension{1});
  name{2}=char(nc.latitude.dimension{2});
  name=sort(name);
  xname=char(name{1});
  yname=char(name{2});
  nc.x=nc.(xname);
  nc.y=nc.(yname);
 end
elseif length(nc.longitude.dimension)>1 || ~strcmp('longitude',nc.longitude.dimension)
 disp(['The longitude dimension is: ',char(nc.longitude.dimension)]);
 error(['Longitude is already 2D: ',char(nc.longitude.dimension)]);
elseif length(nc.latitude.dimension)>1 || ~strcmp('latitude',nc.latitude.dimension)
 disp(['The latitude dimension is:',char(nc.latitude.dimension)]);
 error(['Latitude is already 2D: ',char(nc.latitude.dimension)]);
elseif isfield(nc,'x') | isfield(nc,'y')
 error('x and/or y already variables in the netcdf structure')
end

% add new x and y variables if needed
if calc 

 if debug; disp('Generating generalized grid from a lat-lon grid'); end
 
 % create x and y
 nc=ridgepack_add(nc,'x',[1:length(nc.longitude.data)],'x-coordinate in Cartesian system',...
          {'x'},'','NC_FLOAT');
 nc=ridgepack_add(nc,'y',[1:length(nc.latitude.data)],'y-coordinate in Cartesian system',...
          {'y'},'','NC_FLOAT');

 % Create new variables for latitude and longitude
 lon=nc.longitude.data;
 lat=nc.latitude.data;

 clear nc.longitude.data nc.latitude.data

 [nc.latitude.data,nc.longitude.data]=meshgrat(lat,lon);

 nc.longitude.dimension={'y','x'};
 nc.latitude.dimension={'y','x'};

 xname='longitude';
 yname='latitude';

elseif debug; 
 disp('Reassigning dimension names to 2D lat-lon matrix'); 
end

% Now change dimensions of all other variables to reflect previous change
[nc,variablenames,numbervariables]=ridgepack_sort(nc);
for m = 1:numbervariables
 name=char(variablenames(m));
 for i=1:length(nc.(name).dimension)
  if strcmp(xname,char(nc.(name).dimension{i}))
   nc.(name).dimension{i}='x';
  end
  if strcmp(yname,char(nc.(name).dimension{i}))
   nc.(name).dimension{i}='y';
  end
 end
 if any(strcmp('x',nc.(name).dimension)) && ...
    any(strcmp('y',nc.(name).dimension)) && ...
    ~strcmp(yname,name) && ...
    ~strcmp(xname,name)
  nc.(name).coordinates='latitude longitude';
 end
end


% remove xname and yname variables 
if ~calc
 nc=rmfield(nc,xname);
 nc=rmfield(nc,yname);
end

% sort the structure
nc=ridgepack_sort(nc);

if debug; disp(['...Leaving ',mfilename]); end

