function ridgepack_variable(ncid,nc,rec,recp,recl)

% ridgepack_variable - Adds name, type and attributes of new variable to netcdf file
%
% function ridgepack_variable(ncid,nc,rec,recp)
%
% This function writes a new variable to a netcdf file in CF-1.4 
% convention given:
%   
% Mandatory inputs:
% ncid    - Netcdf ID of the open netcdf file in the working directory.
% nc      - Data structure to be written to netcdf file (see ridgepack_struct)
%
% Optional inputs:
% rec     - Cell array of dimensions for which grid point values are provided
%           in recp if adding a hyperslab to a netcdf file.
%           Here is an example of what rec should look like, if one were 
%           to write a slice of data for a  variable with dimensions 
%           {'time','x','y','z'} on the y and z plane, one would need to 
%           specifiy the values of time and x for which the hyperslab is defined, 
%           hence rec={'time','x'}. The exact indices of time and x that 
%           define the hyperslab must then be provided in recp. 
%           For example recp={3,10} would specify for the 3rd and 10th 
%           indices for time and x, respectively. The unlimited
%           dimension in netcdf files written is reserved only for 
%           a variable named 'time'.  No other variable can take the
%           unlimited dimension using this function.  If time is not
%           a dimension, then an unlimited dimension is not assigned.
%           If you wish to append a time-constant hyperslab to a netcdf
%           file, the you simply specify zero as the recp value for time:
%           ridgepack_write(nc,ncfile,{'time'},{0}).
%
% recp    - Cell array of grid points corresponding to the dimensions in rec
%           for which a hyperslab to be slotted into the netcdf archive.
%
% recl    - Full arrays of non-time dimensions to be written when the netcdf 
%           file is first established. This is required if one intends to append
%           hyperslabs to a file that are not all for a single time slice.
%           If, for example, a variable has dimensions {'x','y','time'}, and one
%           wants to establish a netcdf file, and then write future hyperslabs
%           to that file with for hyperslabs of {'x','y'} (i.e. timeseries
%           for each (x,y) grid point), then would establish the full x and y
%           diemensions in the first write:
%
%           ridgepack_write(ncnew,filename,{'x','y'},{i,j},{xdata,ydata})
%
%           where ncnew only includes data for (x,y)=(i,j), but xdata and ydata
%           includes all values for those dimensions. This establishes the
%           netcdf file. Subsequent writes can be placed inside a loop too append
%           hyperslabs to the file:
%
%           for i=1:length(xdata);
%           for j=1:length(ydata)
%
%            ridgepack_write(ncnew,filename,{'x','y'},{i,j});
%
%           end
%           end
%
%           So in the above example, recl={xdata,ydata}
%
%
% The nc data structure should have the following form for every variable called
% 'name1', 'name2',... ,'nameN' etc for random names to be written to the 
% netcdf file:
%
% nc.nameN.long_name     - description of data
% nc.nameN.units         - units of the data
% nc.nameN.dimension     - cell array with each dimension of 'name'
% nc.nameN.data          - array of values to be written
% nc.nameN.fillvalue     - value assigned to fillvalue and missing_value data 
% nc.nameN.type          - NC_BYTE, NC_CHAR, NC_SHORT, NC_INT, NC_FLOAT or NC_DOUBLE
%
% These are the main fields, but other descriptive fields are also allowed, 
% and these are provided in the documentation for ridgepack_struct ('help ridgepack_struct').
%
% Important points to note in using this function:
%
% 1) Data with several dimensions, for example, would have:
%    nc.nameN.dimension={'time','x','y'} for three dimensions 
%    time, x and y.
%
% 2) If dimension is a single value (not a cell array) and is identical
%    to the name of the variable, then the variable name is established
%    as a dimension in the netcdf file.
%
% 3) data must be an array with an identical number of dimensions as 
%    there are rows in the dimension array and of the same length,
%    sequentially, as the listed dimensions. 
%
% 4) The file ncfile must have already been created with global 
%    attributes using the matlab function ncglobal.
%
% 5) If fillvalue or units are not used or relevant, 
%    they may be omitted. All other fields must be supplied.
%
% 6) If fillvalue is supplied, then it is assumed that all
%    NaNs in the dataset are fill values, and these are substituted
%    with the fillvalue.  If fillvalue is not set and NaNs 
%    exist in the dataset, the function terminates with an error.
%
% 7) The data are not written in this function if they include
%    the time dimension.  In that case, they should be added to the
%    netcdf record using ncnewrecs.m after first defining the 
%    variable using this function.
%
% 8) There is one part of the nc data structure to which these 
%    comments do not apply.  nc.attributes is not used in this
%    function, but is used when initially establishing the file.
%
% 9) When adding a hyperslab to a netcdf file as a timeseries, 
%    the netcdf file and data structure must already be established.
%    See recl above for creating dimensions longer than those initially
%    supplied to the first writing of the file.
%
% The output structure, nc, is altered to include attributes added
% as a result of checking for type and fill values. Also see 
% ridgepack_struct for a description of the netcdf data structure.
%
% Andrew Roberts, Naval Postgraduate School, March 2018  (afrobert@nps.edu)
%

global debug;
if debug; disp(['Entering ',mfilename,'...']); end

% logical to turn on character variable time identifier being added to netcdf file
% to help identify time when using ridgepack_dump.
timestamp=1;

% Run a check on input file name
% check inputs
if nargin==2; 
 hyperslab=0;
 if ~isstruct(nc)
  error('nc input is not a structure')
 elseif ~isnumeric(ncid) | isnan(ncid)
  error('netcdf id is incorrect')
 end
elseif nargin>=4
 hyperslab=1;
 if ~iscell(recp)
  error('recp should be a cell array');
 elseif isempty(recp)
  error('recp entered as an empty array - using slot in the wrong context');
 elseif ~iscell(rec)
  error('rec should be a cell array')
 elseif isempty(rec)
  error('rec entered as an empty array - using slot in the wrong context')
 elseif length(recp)~=length(rec)
  error('recp and rec are different lengths');
 end
 for i=1:length(recp)
  if ~isnumeric(recp{i})
   error(['Element ',num2str(i),' of recp must be numeric'])
  end
  if ~ischar(rec{i})
   error(['Element ',num2str(i),' of rec must be a character string'])
  end
 end
else
 error('Incorrect number of arguments')
end

% set constants
tlen=20;

% Don't prefill the netcdf file 
netcdf.setFill(ncid,'FILL');

% Get number of variables to be written from the input structure.
[nc,variablenames,numbervariables,dimorder]=ridgepack_sort(nc);
xvarnames=ridgepack_cellcat(variablenames);
if length(xvarnames)>50
 disp([num2str(numbervariables),' variables: ',xvarnames(2:50),'... etc']);
else
 disp([num2str(numbervariables),' variables: ',xvarnames(2:end)]);
end

% Get information structure from the netcdf file
[numdims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);

% Set up array to record dimenion variables and unlimited dimensions
diminfo=zeros([1 numbervariables]);

% Check that all necessary dimensions exist and write them if need be
for m = 1:numbervariables

 % Get the name of the variable
 name=char(variablenames(m)); 

 % Check to see if this variable is a dimension and needs to be written.
 if length(nc.(name).dimension)==1 & strcmp(char(nc.(name).dimension{1}),name);

  % record this variable as a dimension,
  if strcmpi('time',name);
   diminfo(m)=1;
  else
   diminfo(m)=2;
  end
  
  % Check that the dimension does not already exist
  writedim=1;
  for i=0:numdims-1
   [dimname,dimlen]=netcdf.inqDim(ncid,i);
   if strcmp(name,dimname);
    if debug; disp([name,' is a dimension.']); end
    writedim=0;
    if dimlen~=length(nc.(name).data) & i~=unlimdimid & ~hyperslab
     dimlen
     length(nc.(name).data)
     error([name,' dimension size incorrect.']);
    end
   end
  end
 
  % Define dimensions, setting time as an unlimited dimension
  if writedim
   netcdf.reDef(ncid);
   if strcmpi('time',name);
    if debug;
     disp(['Setting ',name,' as a dimension with unlimited length']);
     end
    dimid=netcdf.defDim(ncid,'time',netcdf.getConstant('NC_UNLIMITED'));
   else
    % make allowance for hyperslab setup to append timeseries
    if nargin==5 & any(strcmpi(name,rec))
     lengthdim=max(length(nc.(name).data),length(recl{find(strcmpi(name,rec))}));
     disp(['Setting ',name,' as a dimension with length ',num2str(lengthdim)]); 
    else
     lengthdim=length(nc.(name).data);
     if debug; 
      disp(['Setting ',name,' as a dimension with length ',num2str(lengthdim)]); 
     end
    end
    dimid=netcdf.defDim(ncid,name,lengthdim);
   end
   netcdf.endDef(ncid);
  end

 end

end

% Update information structure from the netcdf file
[numdims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
try
 [unlimdimname,unlimdimlen]=netcdf.inqDim(ncid,unlimdimid);
catch
 disp('Unable to find an unlimited dimension')
end

% Now check which variables need to be written
for m = 1:numbervariables

 % Get variable name
 name=char(variablenames(m));

 % Check dimensions of this variable exist in the netcdf file 
 timeexists=0;
 dimid=[];
 if isfield(nc.(name),'dimension')
  for i=1:length(nc.(name).dimension)
   dimname=char(nc.(name).dimension{i});
   try
    dimid(i)=netcdf.inqDimID(ncid,dimname);
   catch
    error(['Unable to find the dimension ',dimname])
   end
   [dimname,dimlen]=netcdf.inqDim(ncid,dimid(i));
   if ~strcmp('time',dimname) & size(nc.(name).data,i) ~= dimlen & ...
      length(nc.(name).dimension)>1 & ~hyperslab
    error(['''',name,''' length for ''',dimname,...
           ''' inconsistent with dimension ''',dimname,'''']);
   end
   
   % check whether time is a dimension, and if so, check for tchar dimension
   if strcmpi(nc.(name).dimension(i),'time') ;
    nc.(name).dimension{i}='time';
    timeexists=i;
   end
  end
 end

 % Check to see if this variable needs to be established in the netcdf file.
 % Do not write dimension variables that are on a set of natural number 
 % indices [1 2 3 4 5 ...]

 try 

  tcharvarid=[];

  % find variables that should only be a dimension and don't check they are a variable
  if ~(length(nc.(name).dimension)==1 && ...
       strcmp(name,char(nc.(name).dimension{1})) && ...
       strcmp('NC_INT',nc.(name).type) && ...
       isfield(nc.(name),'data') && ...
       all(reshape([1:length(nc.(name).data)],size(nc.(name).data))==nc.(name).data))

   varid=netcdf.inqVarID(ncid,name);

   % get the time string variable id
   if ~isempty(findstr(name,'time')) & timestamp & ~strcmpi(name,'time_bounds');
    tcharvarid=netcdf.inqVarID(ncid,[name,'_string']);
   end

  end

  dimwrite=0;

 catch

  % open definitions
  netcdf.reDef(ncid);

  % Set up the new netcdf variable and its attributes from the nc structure 
  if debug; 
   if isempty(dimid)
    disp(['Variable ''',name,''' has no dimension identifiers']);
   else
    disp(['Variable ''',name,''', Dimension identifier(s): ',num2str(dimid)]);
   end
  end


  varid=netcdf.defVar(ncid,name,netcdf.getConstant(nc.(name).type),dimid);
  gatt=fieldnames(nc.(name));
  for i=1:length(gatt)
   attname=char(gatt{i});
   if ~strcmpi(attname,'data') & ...
      ~strcmpi(attname,'weight') & ...
      ~strcmpi(attname,'dimension') & ...
      ~strcmpi(attname,'fillvalue') & ...
      ~strcmpi(attname,'flag_values') & ...
      ~strcmpi(attname,'flag_meanings') & ...
      ~strcmpi(attname,'type') 
    netcdf.putAtt(ncid,varid,attname,nc.(name).(attname));
   end
  end

  % If writing a time variable, add character variable and add a calendar attribute 
  % to time (these are not carried with the nc structure)
  tcharvarid=[];
  if ~isempty(findstr(name,'time'));

   if length(dimid)==1 & ~strcmpi(nc.(name).type,'NC_CHAR') & timestamp

    % create character dimension for time stamp if it does not already exist
    try
     tchardimid=netcdf.inqDimID(ncid,'tchar_len');
    catch
     if debug; disp(['Setting tchar_len as a dimension with length ',num2str(tlen)]); end
     tchardimid=netcdf.defDim(ncid,'tchar_len',tlen);
    end
    [tdimname,tdimlen]=netcdf.inqDim(ncid,tchardimid);
    if tdimlen~=tlen 
     error(['tchar_len is being used by an existing variable with length ',num2str(tlen)])
    end

    % set up the time string variable for name or overwrite existing variable
    try
     tcharvarid=netcdf.defVar(ncid,[name,'_string'],netcdf.getConstant('NC_CHAR'),...
                [tchardimid dimid]);
     netcdf.putAtt(ncid,tcharvarid,'long_name',[name,' values as a string']);
     netcdf.putAtt(ncid,tcharvarid,'units','dd-mmm-yyyy HH:MM:SS');
    catch
     error([name,'_string already exists in the netcdf file']);
    end

   end

   % add the time axis
   if length(nc.(name).dimension)==1 & strcmp(char(nc.(name).dimension),'time');
    netcdf.putAtt(ncid,varid,'axis','T');
   end

  end

  % close definitions
  netcdf.endDef(ncid);

  % If writing a dimension for the first time, the full variable needs to be written
  dimwrite=1;

 end
 
 % Write data for variables with a dimension
 if isempty(nc.(name).dimension);

   if length(nc.(name).data)==1
    netcdf.putVar(ncid,varid,nc.(name).data);
    disp([name,' is dimensionless']) 
   else
    error('Cannot write netcdf file without dimensions for data with length > 1')
   end

 else

   if strcmpi(name,'time') | strcmpi(name,'time_bounds')

    % create time stamp string (first commented line is for adding index)
    if ~isempty(tcharvarid); time_stamp=datestr(nc.(name).data,0); end

    % Check if calendar is already in netcdf file, and adjust data if so
    if  ~isfield(nc.(name),'calendar');
     error('Missing calendar descriptor in nc structure')
    else
     try; 
      calendar=netcdf.getAtt(ncid,varid,'calendar');
      if ~strcmpi(calendar,nc.(name).calendar)
       disp(['Adjusting time data to the ',calendar,' calendar'])   
       nc.(name).calendar=calendar; 
      end
     end
    end

    % Check if time units are already in netcdf file, and adjust data if so
    if  ~isfield(nc.(name),'units');
     error('Missing time units descriptor in nc structure')
    else
     try; 
      units=netcdf.getAtt(ncid,varid,'units');
      if ~strcmpi(units,nc.(name).units)
       disp(['Adjusting time data to ',units])   
       nc.(name).units=units; 
      end
     end
    end

    % Convert time to the chosen units in the ridgepack_structure
    timeshape=size(nc.(name).data);
    [nc.(name).data,units,calen]=ridgepack_timeconvert(nc.(name).data(:),nc.(name).units,...
                                               nc.(name).calendar,1);
    nc.(name).data=reshape(nc.(name).data,timeshape);


    if ~strcmpi(calen,nc.(name).calendar) | ~strcmpi(units,nc.(name).units)
     netcdf.reDef(ncid);
     netcdf.putAtt(ncid,varid,'calendar',calen);
     netcdf.putAtt(ncid,varid,'units',units);
     netcdf.endDef(ncid);
     nc.(name).units=units;
     nc.(name).calendar=calen;
    end

    % Check that times are positive
    if any(nc.(name).data(:)<0)
     disp(['WARNING: Negative time values being written using ',calendar,' ',units])
    end

   end

   % Determine indices of hyperslab and unlimited dimension
   if diminfo(m)<2 & hyperslab

    % get variable ID
    varid=netcdf.inqVarID(ncid,name);

    % get variable information
    [varname,xtype,dimid,natts]=netcdf.inqVar(ncid,varid);
    numdims=length(dimid);

    start=zeros(1, numdims);
    count=zeros(1, numdims);

    if isfield(nc.(name),'dimension') & length(nc.(name).dimension)==1
     varsize=length(nc.(name).data);
    else
     varsize=zeros(size(nc.(name).dimension));
     for ctd=1:length(nc.(name).dimension)
      varsize(ctd)=length(nc.(char(nc.(name).dimension{ctd})).data);
     end
    end

    if debug; disp([name,' varsize = ',num2str(varsize)]); end

    % set start and count vectors against dimension names and lengths
    for i=1:numdims

     [dimname,dimlength] = netcdf.inqDim(ncid,dimid(i));
     ncdimindex=find(strcmp(nc.(name).dimension,dimname));

     if isempty(ncdimindex)
      error(['Unable to find dimension ',dimname,' in the provided nc structure']);
     end

     dimindex=find(strcmp(rec,dimname));

     % this section finds the size of the hyperslab for dimensions not specified in rec
     if isempty(dimindex) 

      % if time is not specified, make the dimension the length of the time data written
      if strcmp(unlimdimname,dimname)
       start(i)=0;
       count(i)=varsize(i);
      else
       start(i)=0;
       count(i)=max(1,dimlength);
      end

      if count(i)~=size(nc.(name).data,ncdimindex) & ndims(nc.(name).data)>=ncdimindex 
       disp([' '])
       disp(['Dimensions disagree between the nc structure and netcdf file:']);
       disp(['-netcdf file ',name,' dimension length: ',num2str(dimlength)])
       disp(['-nc structure ',dimname,' dimension length: ',num2str(size(nc.(name).data,ncdimindex))])
       disp([' '])
       error(['Dimension lengths disagree']);
      end

     % this section appends a new time or unlimited dimension record or records
     elseif dimid(i)==unlimdimid & recp{dimindex}==0

      start(i)=unlimdimlen;
      count(i)=varsize(i);

     % this section simply sets the indices for one slice to be inserted
     else

      start(i)=recp{dimindex}-1;
      count(i)=1;

     end

     if varsize(i)~=count(i)
      varsize
      count
      error(['Size of ',dimname,' in nc structure is inconsistent with the requested hyperslab.'])
     end

    end

   elseif dimwrite | ~hyperslab


    if debug; disp('In: Code segment ''dimwrite | ~hyperslab'''); end

    % write entire data to nc structure
    count=zeros(size(nc.(name).dimension));
    for ctd=1:length(nc.(name).dimension)
     count(ctd)=length(nc.(char(nc.(name).dimension{ctd})).data);
    end
    start=zeros(size(count));

   end

   if dimwrite | diminfo(m)<2

    % get matlab class of variable from equivalent netcdf class
    mclass=ridgepack_onvert(nc.(name).type);

    % check for scale factor and offset and adjust data accordingly
    try
     scale_factor=netcdf.getAtt(ncid,varid,'scale_factor');
    catch
     scale_factor=1;
    end

    try
     add_offset=netcdf.getAtt(ncid,varid,'add_offset');
    catch
     add_offset=0; 
    end

    % assign data to be written, substituting when recl is provided as a fifth argument 
    % to set up adding in hypslabs for more than one time (i.e. along the time axis)
    eval(['output_data=',mclass,'((nc.',name,'.data - add_offset)/scale_factor);']);
    if nargin==5 & any(strcmpi(name,rec))
     if isnumeric(recl{find(strcmpi(name,rec))}) | ischar(recl{find(strcmpi(name,rec))})
      if count<length(recl{find(strcmpi(name,rec))})
       eval(['output_data=',mclass,'((recl{find(strcmpi(''',name,...
             ''',rec))} - add_offset)/scale_factor);'])
       count=length(recl{find(strcmpi(name,rec))})
      end
     else
      error(['recl element matching ',name,' is neither a numeric nor a character array']) 
     end
    end

    % Scale, offset and convert data to the type to be written
    if ~all(isnan(nc.(name).data(:)))
     if isfield(nc.(name),'valid_range') 
      if ~isnumeric(nc.(name).valid_range) 
       error(['valid range for ',name,' should be a numeric value [min max]'])
      elseif length(nc.(name).valid_range)~=2 
       error(['valid range for ',name,' have length 2 with value [min max]'])
      else
       valrange=sort(nc.(name).valid_range);
      end
     else
      valrange=[min(nc.(name).data(~isnan(nc.(name).data))) ...
 	        max(nc.(name).data(~isnan(nc.(name).data)))];
     end
     eval(['valrange=',mclass,'(valrange);']);
     novalrange=false;
    else
     novalrange=true;
    end

    % check for fill values and fill the data where needed
    if isfield(nc.(name),'fillvalue') 
      netcdf.reDef(ncid);
      if ~novalrange; netcdf.putAtt(ncid,varid,'valid_range',valrange); end
      eval(['fillval=',mclass,'(nc.(name).fillvalue);']); 
      netcdf.putAtt(ncid,varid,'_fillvalue',fillval);
      netcdf.endDef(ncid);
      output_data(isnan(nc.(name).data))=fillval;
    elseif strcmp(nc.(name).type,'NC_INT') | strcmp(nc.(name).type,'NC_BYTE')
      netcdf.reDef(ncid);
      if ~novalrange; netcdf.putAtt(ncid,varid,'valid_range',valrange); end
      netcdf.endDef(ncid);
    end

    % check for flag values
    if isfield(nc.(name),'flag_values') & isfield(nc.(name),'flag_meanings')
      netcdf.reDef(ncid);
      eval(['nc.',name,'.flag_values=',mclass,'(nc.',name,'.flag_values);']);
      netcdf.putAtt(ncid,varid,'flag_values',nc.(name).flag_values);
      netcdf.putAtt(ncid,varid,'flag_meanings',nc.(name).flag_meanings);
      netcdf.endDef(ncid);
    elseif ~isfield(nc.(name),'flag_values') & isfield(nc.(name),'flag_meanings')
      error('missing flag values')
    elseif isfield(nc.(name),'flag_values') & ~isfield(nc.(name),'flag_meanings')
      error('missing flag meanings')
    end

    % write the value
    if debug;
     disp(['Writing ',name,': start=[',num2str(start),'] count=[',num2str(count),']']);
     disp(['Size of data placed into netcdf = ',num2str(size(output_data))]); 
    end
    netcdf.putVar(ncid,varid,start,count,output_data);

    % write the auxiliary time stamp string if needed
    if ~isempty(tcharvarid); 
      start=[0 start];
      count=[tlen count];
      if debug; disp(['Writing ',name,'_string: start=[',num2str(start),...
                      '] count=[',num2str(count),']']); end
      netcdf.putVar(ncid,tcharvarid,start,count,time_stamp'); 
    end

   end

 end

end % for m = 1:numbervariables

if debug; disp(['...Leaving ',mfilename]); end


