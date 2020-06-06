function [nc]=ridgepack_coords(nc,name);

% ridgepack_coords - Assigns coordinates attribute within a netcdf structure
%
% function [nc]=ridgepack_coords(nc,name);
%
% This function assigns the coordinates attribute to a netcdf structure if it is 
% appropriate. It also changes the names of lat/long variables to latitude and 
% longitude, where appropriate, based on the coordinates found in the netcdf structure. 
%
% INPUT:
%
% nc - netcdf structure (see ridgepack_struct for more information)
% name - optional input to assign the coordinates for a specific variable
%
%
% OUTPUT:
%
% nc - netcdf structure with latitude and longitude variable names, and coordinate
%      attributes where a variable has the same dimension variable as latitude and 
%      longitude.
%
% Ridgepack Version 1.0
% Andrew Roberts, Naval Postgraduate School, March 2018 (afrobert@nps.edu)
% Naval Postgraduate School, Department of Oceanography

global debug;
if debug; disp(['Entering ',mfilename,'...']); end

% set up function inputs if not specified
if nargin<1 | nargin>2
 error('ridgepack_coords: Must have at one nc structure input and an option variable name')
elseif nargin>=1 & ~isstruct(nc); 
 error('ridgepack_coords: nc is not a structure')
end

% quick sort
[nc,variablenames,numbervariables]=ridgepack_sort(nc);

% assign coordinates to the entire nc structure
[dimc,coor,dims,vars,mdim,mesh]=ridgepack_content(nc);


% Locate variable from input argument if provided
if nargin==2 & ~isempty(name)
 if ischar(name)
  if any(strcmp(name,vars)) & isfield(nc.(name),'coordinates')
   varc=name;
  else
   varc=[];
  end
 else
  error('ERROR: Name must be a character variable')
 end
else
 varc=[];
end

% Otherwise check to see if data variable coordinates 
% are identical throughout the vars list, and hence there
% is only one possible coordinate descriptor choice.
coordv=[];
if isempty(varc) & length(vars)>=1
 for m = 1:length(vars)
  varc=char(vars{m});
  if isfield(nc.(varc),'coordinates') 
   if isempty(coordv)
    coordv=nc.(varc).coordinates;
   elseif ~strcmp(coordv,nc.(varc).coordinates)
    coordv=['n/a'];
   end
  end
 end
 if strcmp(coordv,'n/a')
  varc=[];
  coordv=[];
 end
end


% Process these options if varc not found
if any(strcmp('latitude',variablenames)) & ...
   any(strcmp('longitude',variablenames))

  found=false;

  map_def=[];
  for m = 1:numbervariables
    name=char(variablenames(m));
    if isfield(nc.(name),'grid_mapping_name')
     map_def=name;
    end
  end

  for m = 1:numbervariables
   name=char(variablenames(m));
   if any(strcmp('latitude',nc.(name).dimension)) & ...
      any(strcmp('longitude',nc.(name).dimension)) 
     if ~isfield(nc.(name),'coordinates')
      nc.(name).coordinates='latitude longitude';
      if debug; disp(['Adding coordinates attribute to ',name]); end
      found=true;
     end
     if ~isempty(map_def) & ~isfield(nc.(name),'grid_mapping') 
      nc.(name).grid_mapping=map_def;
      if debug; disp(['Adding grid mapping attribute to ',name]); end
     end
   end
  end

  if found; return; end

  if isempty(setdiff(nc.latitude.dimension,nc.longitude.dimension))
   for m = 1:numbervariables
    name=char(variablenames(m));
    if ~strcmp('latitude',name) & ...
       ~strcmp('longitude',name) & ...
        length(nc.(name).dimension)>1 & ...
        isempty(intersect({'time'},nc.latitude.dimension)) & ...
        isempty(intersect({'time'},nc.longitude.dimension)) & ...
       ~isfield(nc.(name),'coordinates') & ...
       ~any(strcmpi(name,coor))

	if isempty(setdiff(union(nc.latitude.dimension, ...
                   nc.(name).dimension),nc.latitude.dimension)) | ...
	   strcmpi(setdiff(union(nc.latitude.dimension, ...
                   nc.(name).dimension),nc.latitude.dimension),'time')

	  	nc.(name).coordinates='latitude longitude';
	  	if debug; disp(['Assigned new coordinates attribute for ',name]); end

		if ~isempty(map_def) 
		 nc.(name).grid_mapping=map_def;
	  	 disp(['Assigned new mapping attribute for ',name])
		end

	end
    end
   end
  end

  return

end

% find all coordinate variables and possible lat/lon variables
allcoords={};
lonname={};
latname={};
lonlongname={};
latlongname={};
for m = 1:numbervariables
  name=char(variablenames(m));
  coords=ridgepack_varcoords(nc,name);
  allcoords=union(coords,allcoords);
  if (~isempty(strfind(name,'lon')) | ~isempty(strfind(name,'Lon')) | ~isempty(strfind(name,'LON'))) & isempty(strfind(name,'longwave')) & isempty(strfind(name,'Longwave'))
   lonname=union(lonname,{name});
   if ~isempty(strfind(nc.(name).long_name,'longitude')) | ~isempty(strfind(nc.(name).long_name,'Longitude'))
    lonlongname=union(lonlongname,{name});
   end
  end 
  if ~isempty(strfind(name,'lat')) | ~isempty(strfind(name,'Lat')) | ~isempty(strfind(name,'LAT')) 
   latname=union(latname,{name});
   if ~isempty(strfind(nc.(name).long_name,'latitude')) | ~isempty(strfind(nc.(name).long_name,'Latitude'))
    latlongname=union(latlongname,{name});
   end
  end
end


newlonname=intersect(lonlongname,lonname);
newlatname=intersect(latlongname,latname);
if ~isempty(newlatname); latname=newlatname; end
if ~isempty(newlonname); lonname=newlonname; end


% account for no latitudes or longitudes
if isempty(lonname) | isempty(latname)

  if debug; disp('No 2D latitude/longitude variables found in this nc structure'); end
  return

% assign latitude and longitude references if there is no dispute
elseif length(lonname)==1 & length(latname)==1

  oldlon=char(lonname);
  oldlat=char(latname);

  disp(['Latitude and longitude variables reassigned from ',oldlat,' and ',oldlon]);

% else search for latitude and longitude references if in coords
elseif length(union(allcoords,lonname))==1 & ...
       length(union(allcoords,latname))==1

  oldlon=char(union(allcoords,lonname));
  oldlat=char(union(allcoords,latname));

  disp(['Latitude and longitude coordinates reassigned from ',oldlat,' and ',oldlon]);

% else search for latitude and longitude references if in coordinates attribute
elseif ~isempty(varc) && isfield(nc.(varc),'coordinates')

  for i=1:length(latname)
   if ~isempty(strfind(nc.(varc).coordinates,char(latname{i})))
    oldlat=char(latname{i});  
   end
  end

  for i=1:length(lonname)
   if ~isempty(strfind(nc.(varc).coordinates,char(lonname{i})))
    oldlon=char(lonname{i});  
   end
  end

  disp(['Latitude and longitude coordinates reassigned from ',oldlat,' and ',oldlon]);

% if there are multiple options
else

 for i=1:length(lonname)
	 disp([num2str(i),': ',char(lonname{i})]);
 end
 lonc=0;
 while ~(isnumeric(lonc) & (lonc>=1 & lonc<=length(lonname)))
  try
   lonc=input('Enter longitude coordinate (not dimension) index: ');
   if ~(lonc>=1 & lonc<=length(lonname)) | ~isnumeric(lonc) ; error ; end
  catch
   lonc=0;
   disp('Incorrect input')
  end
 end
 oldlon=char(lonname{lonc});

 for i=1:length(latname)
	 disp([num2str(i),': ',char(latname{i})]);
 end
 latc=0;
 while ~(isnumeric(latc) & (latc>=1 & latc<=length(latname)))
  try
   latc=input('Enter latitude coordinate (not dimension) index: ');
   if ~(latc>=1 & latc<=length(latname)) | ~isnumeric(latc) ; error ; end
  catch
   latc=0;
   disp('Incorrect input')
  end
 end
 oldlat=char(latname{latc});

end

% reassign variable name, changing dimensions and coordinate names and 
% looking for a map definition
map_def=[];
for m = 1:numbervariables
  name=char(variablenames(m));

  % replace coordinates descriptor
  if isfield(nc.(name),'coordinates')
   nc.(name).coordinates=regexprep(char(nc.(name).coordinates),oldlon,'longitude');
   nc.(name).coordinates=regexprep(char(nc.(name).coordinates),oldlat,'latitude');
  end

  if isfield(nc.(name),'grid_mapping_name')
   map_def=name;
  end

  if strcmp(name,oldlon) & ~strcmp('longitude',name)
   nc.longitude=nc.(name);
   nc=rmfield(nc,name);
   name='longitude';
   variablenames{m}='longitude';
  end

  for i=1:length(nc.(name).dimension)
   if strcmp(oldlon,char(nc.(name).dimension{i}))
	   nc.(name).dimension{i}='longitude';
   end
  end

  if strcmp(name,oldlat) & ~strcmp('latitude',name)
   nc.latitude=nc.(name);
   nc=rmfield(nc,name);
   name='latitude';
   variablenames{m}='latitude';
  end

  for i=1:length(nc.(name).dimension)
   if strcmp(oldlat,char(nc.(name).dimension{i}))
	   nc.(name).dimension{i}='latitude';
   end
  end

end

% Correctly assign coordinate variables
if isfield(nc,'latitude') & isfield(nc,'longitude')
 if isempty(setdiff(nc.latitude.dimension,nc.longitude.dimension))
   for m = 1:numbervariables
    name=char(variablenames(m));
    if ~strcmp('latitude',name) & ...
       ~strcmp('longitude',name) & ...
        length(nc.(name).dimension)>1 & ...
        isempty(intersect({'time'},nc.latitude.dimension)) & ...
        isempty(intersect({'time'},nc.longitude.dimension)) & ...
        ((~isfield(nc.(name),'coordinates') & ~any(strcmpi(name,coor))) | ...  
        (isfield(nc.(name),'coordinates') &&  strcmp(nc.(name).coordinates,[oldlat,' ',oldlon])))

	if isempty(setdiff(union(nc.latitude.dimension, ...
           nc.(name).dimension),nc.latitude.dimension)) | ...
	   strcmpi(setdiff(union(nc.latitude.dimension, ...
           nc.(name).dimension),nc.latitude.dimension),'time')
      
	  	nc.(name).coordinates='latitude longitude';
	  	if debug; disp(['Assigned new coordinates attribute for ',name]); end

		if isfield(nc,'cell_area')
			nc.(name).cell_measures='area: cell_area';
		end

		if ~isempty(map_def) & ~isfield(nc.(name),'grid_mapping')
			nc.(name).grid_mapping=map_def;
		end

	end
    end
   end
 end
end


% Remove time from the coordinates if it has been put there
removet=0;
for m = 1:numbervariables
  name=char(variablenames(m));
  coords=setdiff(ridgepack_varcoords(nc,name),{'time'});
  if ~isempty(coords)
   nc.(name).coordinates=ridgepack_cellcat(coords);
   removet=1;
  end
end

if removet & debug; disp(['Removed time from coordinate descriptors.']); end

if debug; disp(['...Leaving ',mfilename]); end


