function [nc,out]=ncstruct(nc,varname)

% NCSTRUCT - Groom, sort and summarize a netcdf structure
%
% function [nc,out]=ncstruct(nc,varname)
%
% This function checks, grooms and sorts the netcdf structure
% nc, then lists the content of the structure, including 
% basic statistics about each netcdf variable. If varname
% is provided (optional), only details about this variable are 
% provided, along with all values from that variable.
%
% Input:
% nc      - netcdf structure (see explanation below).
% varname - variable name from the structure (optional)
%
% Output:
% nc      - groomed and sorted netcdf structure.
% out     - this argument is included only if summary output
%           is required to be written to screen, rather than
%           all information.
%
% Explanation of the netcdf structure
% -----------------------------------
% In this explanation, 'netcdf structure' and 'nc structure' are equivalent.
%
% nc is a structure with the following components:
%
% nc.attribute - global attributes
% nc.name1     - variable object
% nc.name2     - second variable object
% ...
% nc.nameN     - Nth variable object
%
% This can be easily written to or read from a netcdf file
% and is a useful way of storing and handling data within
% matlab, along with the timeseries structure. The global
% attributes of the structure have the following components:
%
% nc.attributes.title - A succinct description of what is in the dataset.
%                       This is compulsory.
%
% nc.attributes.source - The method of production of the original data. 
%           If it is model-generated, source should name the model and its,
%           version as specifically as is useful. If it is an observational, 
%           source should characterize this (e.g., "surface observation" or 
%           "radiosonde"). This attribute should also include specifics of
%           of the instrumentation used to produce the data:
%           1) satellite/remotely sensed data: include summary 
%              of the satellite, instrument and processing algorithms.
%           2) GPS/positional data: include processing procedure and 
%              retrieval method.
%           3) directly sensed data: exact details on sensors used (model
%              numbers, sensor types)
%           4) model generated data: computer on which the data was 
%              generated.
%           5) conglomerate dataset (such as climatology or assimilated):
%              list (within reasonable bounds) all sources of data and 
%              methods/models used to compile the dataset.
%           This is compulsory.
%
% nc.attributes.references - Published or web-based references that describe 
%           the data or methods used to produce it. This need to be included 
%           if references are not applicable.  Exclusion of this will
%           will remove the entry from the netcdf file. This is optional.
%
% nc.attributes.comment - Miscellaneous information about the data or
%           methods use to produce it.  This is optional.
%
% If attribute text is long, it is advisable to place newline markers ("\n") 
% at appropriate places within the text. Each variable component of the netcdf 
% structure normally contains the following for the Nth variable name:
% 
% nc.nameN.long_name - text string of the long version of the variable name
%
% nc.nameN.units  - text string specifying the data units.  If name='time'
%                   then the text string represents the conversion 
%                   to the CF coordinate convention from matlan serial
%                   time when the time dimension data is written to netcdf.
%                   If no time units are specified, the default is
%                   hours since 1-1-1900 00:00:00.0 UTC.
%
% nc.nameN.dimension - cell arrays of texts specifying the dimensions of the 
%                   the variables. So the variables time would have
%                   nc.name.dimension={'time'}, as would a time series of, 
%                   for example, barometric pressure. Multidimensional
%                   data would have multiple dimensions, such as MSL from
%                   the ERA40 dataset: {'time','latitude','longitude'}
%                   where time, latitude and longitude are also defined
%                   variables in the netcdf structure.
%                          
% nc.nameN.data   - data values with the dimension allocated in line with 
%                   the defined dimensions.
%
% nc.nameN.fillvalue - the fill value of the data to be assigned during 
%                   write of the data to netcdf.  The fill value
%                   is assumed to be NaN in matlab, and so all NaN's are
%                   substituted for nc.name.fillvalue when writing to netcdf.
%
% nc.nameN.type   - type the data is to be assigned to when writing to netcdf.
%                   The data is stored in the structure nominally as double
%                   precision, but this can be reduced to single upon writing
%                   to netcdf.
%
% nc.nameN.coordinates - attribute describing which variables define spatial 
%                   coordinates for this variable.
%
% nc.nameN.grid_mapping - convention depicting the type of mapping procedure
%                   such polar stereographic or spherical coordinates.
%
% otherwise the variable may contain metadata only, such as the grid mapping, 
% in which case the variable will typically of type NC_CHAR and will follow 
% CF conventions, or combinations of the above fields and other CF convention
% fields not mentioned above..
%
% An example of output from ncstruct showing a netcdf structure of ERA40 data
% ---------------------------------------------------------------------------
%
% netcdf structure nc {
%  attributes {global}
%     Conventions: 'CF-1.1'
%         history: '2006-02-28 02:34:30 GMT by mars2netcdf-0.92'
%  
%  4 variables, dimension=[ ], data=( )
%  
%  longitude [0 to 357.5, range 357.5]
%     long_name: 'longitude'
%         units: 'degrees_east'
%     dimension: {'longitude'}
%          data: [144x1 double]
%          type: 'NC_FLOAT'
%  
%  latitude [-90 to 90, range 180]
%     long_name: 'latitude'
%         units: 'degrees_north'
%     dimension: {'latitude'}
%          data: [73x1 double]
%          type: 'NC_FLOAT'
%  
%  time [01-Jan-1996 00:00:00 to 30-Jun-1997 18:00:00, mean steps 06:00:00.000]
%     long_name: 'time'
%         units: 'hours since 1900-01-01 00:00:0.0'
%     dimension: {'time'}
%          data: [2188x1 double]
%          type: 'NC_FLOAT'
%  
%  msl (92941.25 to 107756.1875, range 14814.9375, median 101168.6812)
%     long_name: 'Mean sea level pressure'
%         units: 'Pa'
%     dimension: {'time'  'latitude'  'longitude'}
%          data: [2188x73x144 double]
%     fillvalue: -32767
%          type: 'NC_FLOAT'
% }
%
% Tools for interrogating and manipulating data in nc structures
% --------------------------------------------------------------
% For an explanation of tools for nc structures, see ncpolar
% (type 'help ncpolar').
% 
% Written by Andrew Roberts 2012
% Department of Oceanography, Naval Postgraduate School
%
% $Id$

global debug;
if debug; disp(['Entering ',mfilename,'...']); end

% check for latitude and longitude variables to assign
nc=nccoords(nc);

% sort the structure
[nc,variablenames,numbervariables,dimorder]=ncsort(nc);

% Pass through the structure to make sure dimensions and data size are
% consistent by adding leading singleton dimensions where they are missing.
for m = 1:numbervariables
 name=char(variablenames(m));
 if isfield(nc.(name),'data') && isnumeric(nc.(name).data)
  if length(nc.(name).dimension)>ndims(nc.(name).data)
   for i=1:length(nc.(name).dimension)
    if length(nc.(nc.(name).dimension{i}).data)==1 & ...
      ndims(nc.(name).data)>=i & size(nc.(name).data,i)~=1
       resh=size(nc.(name).data);
       resh(i+1:end+1)=resh(i:end);
       resh(i)=1;
       nc.(name).data=reshape(nc.(name).data,resh);
    end
   end
  end
 end
end


% shuffle dimensions into the the order as listed in nc;
nc=ncshuffle(nc,dimorder);

disp(' ');

if nargout==2; 
	out='Summary output provided to screen' ; 
	disp([inputname(1),' summary {']);
else
	disp(['netcdf structure ',inputname(1),' {']);
end

% preamble for complete data summary only
if nargin==1  & nargout<=1;
 if isfield(nc,'attributes');
  disp(' attributes {global}');
  disp(nc.attributes);
  disp(' ');
 else
  disp('no global attributes.');
  disp(' ');
 end
 disp([' ',num2str(numbervariables),...
       ' variables, Key: descriptor={ }, dimension=[ ], data=( )']);
 disp(' ');
end

if not(strmatch('time',variablenames,'exact'))  & nargout<=1;
 disp('There are no recognized time coordinates in this data');
end

% Output information for each variable
grid_map=false;
foundvar=false;
for m = 1:numbervariables

  name=char(variablenames(m));

  if nargin==1 || strcmp(name,varname);

   foundvar=true;

   if isfield(nc.(name),'long_name') & ...
      isfield(nc.(name),'dimension') & ...
      isfield(nc.(name),'data') & ...
      ~isempty(nc.(name).long_name) & ...
      ~isempty(nc.(name).dimension) & ...
      ~isempty(nc.(name).data) ;

    if strcmp(char(nc.(name).dimension{1}),name) &...
       length(nc.(name).dimension)==1;

     if nargin>1 & nargout<=1; disp('DIMENSION:');disp(' ');end

     if strcmpi(name,'time');

      % make sure time array is 1-D
      if length(nc.(name).data)~=length(nc.(name).data(:));;
	      error('time array must be one-dimensional');
      end

      % filter out silly times exceeding 2300
      nc.(name).data(nc.(name).data>datenum(2300,1,1))=NaN;

      % calculate the time steps in the data if record > 1
      l=length(nc.(name).data);
      if(l>1);
       timesteps=nanmedian(nc.(name).data(2:l)-nc.(name).data(1:l-1));
       if timesteps>1;
        steps=['~',num2str(timesteps),' days'];
       elseif timesteps==1;
        steps=[num2str(timesteps),' day'];
       else
        steps=datestr(timesteps,'HH:MM:SS.FFF');
       end
       disp([' ',name,' [',datestr(min(nc.(name).data(~isnan(nc.(name).data)))),...
        	    ' to ',datestr(max(nc.(name).data(~isnan(nc.(name).data)))),...
                    ', timestep ',steps,', length ',num2str(l),']']);
      else
       disp([' ',name,' [',datestr(min(nc.(name).data(~isnan(nc.(name).data)))),...
                    ', length ',num2str(l),']']);
      end

     else

      l=length(nc.(name).data);
      steps=nanmedian(nc.(name).data(2:l)-nc.(name).data(1:l-1));

      if isnan(steps) | isempty(steps)
       disp([' ',name,' [',num2str(min(nc.(name).data(:))),', length ',num2str(l),']']);
      else
       disp([' ',name,' [',num2str(min(nc.(name).data(:))),...
             ' to ',num2str(max(nc.(name).data(:))),...
             ', range ',num2str(range(nc.(name).data(:))),...
             ', steps ',num2str(steps),', length ',num2str(l),']']);
      end

     end

    else	

     if nargin>1 & nargout<=1; disp('DATA:');disp(' ');end

     if length(nc.(name).dimension)>1

      dimchar1=char(nc.(name).dimension{1});
      dimchar2=char(nc.(name).dimension{2});

      if any(strcmp(name,{'longitude','longitude_corner'})) 
       if strcmp(dimchar1(1),'x') 
	       steps=abs(nanmedian(nc.(name).data(1:end-1,1)-nc.(name).data(2:end,1)));
       elseif strcmp(dimchar2(1),'x')
	       steps=abs(nanmedian(nc.(name).data(1,1:end-1)-nc.(name).data(1,2:end)));
       else
        steps=[];
       end
      elseif any(strcmp(name,{'latitude','latitude_corner'})) 
       if strcmp(dimchar1(1),'y') 
	       steps=abs(mean(nc.(name).data(2:end-2,1)-nc.(name).data(3:end-1,1)));
       elseif strcmp(dimchar2(1),'y')
	       steps=abs(mean(nc.(name).data(1,2:end-2)-nc.(name).data(1,3:end-1)));
       else
        steps=[];
       end
      else
       steps=[];
      end

     else

       steps=[];

     end

     if nargout>1 & isfield(nc.(name),'units') & ~strcmpi(name,'time_bounds'); 
	dispdim=[' ',nc.(name).units,'):',nccellcat(nc.(name).dimension)]; 
     elseif nargout>1 
	dispdim=['):',nccellcat(nc.(name).dimension)]; 
     else; 
	dispdim=')'; 
     end

     if strcmpi(name,'time') | strcmpi(name,'time_bounds');

      % make sure time array is 1-D
      if ~strcmpi(name,'time_bounds') & ...
         length(nc.(name).data)~=length(nc.(name).data(:));;
	      error('time array must be one-dimensional');
      end

      % calculate the time steps in the data if record > 1
      l=length(nc.(name).data);
      if(l>1);
       timesteps=nanmedian(nc.(name).data(2:l)-nc.(name).data(1:l-1));
       if timesteps>1;
        steps=[num2str(timesteps),' days'];
       elseif timesteps==1;
        steps=[num2str(timesteps),' day'];
       else
        steps=datestr(timesteps,'HH:MM:SS.FFF');
       end
       if strcmpi(name,'time_bounds')
        disp([' ',name,' (',datestr(min(nc.(name).data(:)),0),...
        	    ' to ',datestr(max(nc.(name).data(:)),0),dispdim])
       else
        disp([' ',name,' (',datestr(min(nc.(name).data(:)),0),...
        	    ' to ',datestr(max(nc.(name).data(:)),0),...
                    ', step ',steps,', ',dispdim]);
       end
      else
       disp([' ',name,' (',datestr(min(nc.(name).data),0),dispdim]);
      end

     elseif isempty(steps) & nargout>0 & prod(size(nc.(name).data))>16000000
      disp([' ',name,' (statistics withheld for efficiency):',...
	   nccellcat(nc.(name).dimension)]); 

     elseif isempty(steps) | isnan(steps)
      disp([' ',name,' (',num2str(min(nc.(name).data(:))),...
            ' to ',num2str(max(nc.(name).data(:))),...
            ...%', range ',num2str(range(nc.(name).data(:))),...
            ', median ',num2str(nanmedian(nc.(name).data(:))),dispdim]);

     else
      disp([' ',name,' (',num2str(min(nc.(name).data(:))),...
            ' to ',num2str(max(nc.(name).data(:))),...
            ...%', range ',num2str(range(nc.(name).data(:))),...
            ', delta ',num2str(steps),dispdim]);
     end

    end

    if nargout<=1;
	    disp(nc.(name));
	    disp(' ');
    end

    if nargin>1 & nargout<=1;
     if strcmp(name,'time');
      for i=1:3:length(nc.time.data)
       if i==length(nc.time.data)
        disp([num2str(i),': ',datestr(nc.time.data(i))])
       elseif i==length(nc.time.data)-1
        disp([num2str(i),': ',datestr(nc.time.data(i)),',  ',...
              num2str(i+1),': ',datestr(nc.time.data(i+1))]);
       else
        disp([num2str(i),': ',datestr(nc.time.data(i)),',  ',...
              num2str(i+1),': ',datestr(nc.time.data(i+1)),',  ',...
              num2str(i+2),': ',datestr(nc.time.data(i+2))]);
       end
      end
     else
      disp(num2str(nc.(name).data(:)));
     end
    end

   elseif isfield(nc.(name),'dimension') && ...
          isempty(nc.(name).dimension) && ...
          strcmp('NC_CHAR',nc.(name).type) ; 

    if isfield(nc.(name),'grid_mapping_name')
     if grid_map
      error('Already encountered one grid mapping for this structure');
     else
      disp([' ',name,' {grid mapping}']);
      grid_map=true;
     end
    else
     disp([' ',name])
    end

    if nargout<=1;
     disp(rmfield(nc.(name),'dimension'));
     disp(' ');
    end

   else

    if ~isfield(nc.(name),'data') || isempty(nc.(name).data); 
	disp([' ',name,': WARNING: MISSING DATA']); 
        nc.(name).data=[];
    end

    if ~isfield(nc.(name),'dimension') || isempty(nc.(name).dimension); 
	disp([' ',name,': WARNING: MISSING DIMENSION']); 
        nc.(name).dimension={};
    end
    if ~isfield(nc.(name),'long_name') || isempty(nc.(name).long_name); 
	nc.(name).long_name=name;
    	disp(name);
    end
    disp(nc.(name));
    disp(' ');

   end

  end

end

if not(foundvar) ; 
	error(['No variables found named ',varname]);
end

disp('}');
disp(' ');

if nargout==0; clear nc; end

if debug; disp(['...Leaving ',mfilename]); end

