function [nc,variablenames,numbervariables,dimorder]=ridgepack_sort(nc)

% ridgepack_sort - Reorders and grooms a netcdf structure
% 
% function [nc,variablenames,numbervariables,dimorder]=ridgepack_sort(nc)
%
% Input:
% nc - netcdf data structure defined in ridgepack_struct.
%
% Output:
% nc              - Cleaned netcdf data structure
% variablenames   - cell array of all variables names in nc
% numbervariables - number of cells in variablenames
% dimorder        - order of dimensions listed in ridgepack_structure
%
% The netcdf data structure is sorted into the 
% following order:
%
%  {
%  nc.attributes
%  nc.name1 (variable defining a dimension)
%  ..
%  nc.namel (variable defining a dimension)
%  nc.namen (data variable)
%  ...
%  nc.namem (data variable)
%  }
%
% Andrew Roberts, Naval Postgraduate School, March 2018  (afrobert@nps.edu)
%

global debug;
if debug; disp(['Entering ',mfilename,'...']); end

% Check data is a structure
if not(isstruct(nc));
	error([inputname(1),' is not a structure']);
end     

% Get names of variables within the netcdf structure.
[variablenames,numbervariables]=ridgepack_name(nc);

% Groom each field for each variable in the structure.
for m = 1:numbervariables
  name=char(variablenames(m));
  nc=ridgepack_groom(nc, name);
end


% Extract dimension and descriptor information, as well as
% dimension specific to character variables.  Add an empty dimension entry
% for character variables without one (to allow for, for example, grid descriptions)
c={};
de={};
chdim={};
nochdim={};
for m = 1:numbervariables
  name=char(variablenames(m));
  if isfield(nc.(name),'type') && strcmp(nc.(name).type,'NC_CHAR') && ~isfield(nc.(name),'dimension')
   nc.(name).dimension={};
  elseif ~isfield(nc.(name),'dimension')
   nc.(name).dimension={};
  end
  len=length(nc.(name).dimension);
  if len>0
   for n = 1:len
    c=union(c,{char(nc.(name).dimension{n})});
   end
  else
    de=union(de,{name});
  end
  if ~isfield(nc.(name),'type')
   nc.(name).type='NC_CHAR';
   if isfield(nc.(name),'data')
    if strcmpi(class(nc.(name).data),'double')
     nc.(name).type='NC_DOUBLE';
    elseif strcmpi(class(nc.(name).data),'int32')
     nc.(name).type='NC_INT';
    elseif strcmpi(class(nc.(name).data),'int64')
     nc.(name).type='NC_INT';
    elseif strcmpi(class(nc.(name).data),'int8')
     nc.(name).type='NC_BYTE';
    elseif strcmpi(class(nc.(name).data),'uint8')
     nc.(name).type='NC_BYTE';
    elseif strcmpi(class(nc.(name).data),'single')
     nc.(name).type='NC_FLOAT';
    end
   end
  end
  if strcmp(nc.(name).type,'NC_CHAR')
   for n = 1:len
    chdim=union(chdim,{char(nc.(name).dimension{n})});
   end
  elseif ~strcmp(name,nc.(name).dimension)
   for n = 1:len
    nochdim=union(nochdim,{char(nc.(name).dimension{n})});
   end
  end
end

% Find the character length dimensions
chdim=setdiff(chdim,nochdim);

% sort the variables to put dimensions after attributes and descriptors
if isfield(nc,'attributes');
	s.attributes=1;
else
	nc.attributes.title=input('Enter a title for the nc structure: ','s');
	s.attributes=1;
end

% add descriptors
if ~isempty(de)
	for m=1:length(de)
		s.(char(de{m}))=1;
	end
end

% sort dimensions in the order here, with character lengths last, 
% then alphabetically for other dimensions
if isfield(nc.attributes,'source') && ...
   strcmpi(nc.attributes.source,'NCEP/DOE AMIP-II Reanalysis (Reanalysis-2) Model') % NCEP
 dimorder=[{'time','latitude','longitude','z','y','y_corner','x','x_corner',...
            'level','lev','depth','height'},chdim];
else % default
 dimorder=[{'time','z','nc','nj','y','y_corner','ni','x','x_corner',...
            'level','lev','depth','height','latitude','longitude'},chdim];
end

% check for character dimensions, and place them last
for i=1:length(dimorder)
 if ~isempty(find(strcmp(c,char(dimorder{i})))); s.(char(dimorder{i}))=1; end
end
for m = 1:length(c)
 s.(char(c{m}))=1;
end

% sort variable first alphabetically
clear c;
l=0;
for m = 1:numbervariables
	name=char(variablenames(m));
	if not(strcmp(name,'exact'));
		l=l+1;
		c{l}=name;
	end
end
c=sort(c);

% Now sort according to the number of dimensions
for d = 0:6
 for l = 1:length(c);
	if length(nc.(char(c{l})).dimension)==d
		s.(char(c{l}))=1;
        end
	if length(nc.(char(c{l})).dimension)>6
		error(['Too many dimensions for ',char(c{l})]);
	end
 end
end

% Check that the same number of dimension variables were obtained
% as the number of dimensions found
[names,numbers]=ridgepack_name(s);
if numbers ~= numbervariables ; 
	disp(['number of expected variables=',num2str(numbers)]);
	disp(['number of variables=',num2str(numbervariables)]);
	disp(['Expected variables: ',ridgepack_cellcat(names)]);
	disp(['Actual variables: ',ridgepack_cellcat(variablenames)]);
	disp('ERROR: Variables and dimensions are not coherent in the nc structure');
        nc
	error('You may be missing some dimension variables in your nc structure');
end

% reorder the fields based on the s structure
nc=orderfields(nc,s);

% Update names of variables within the newly ordered netcdf structure.
[variablenames,numbervariables]=ridgepack_name(nc);

if debug; disp(['...Leaving ',mfilename]); end

