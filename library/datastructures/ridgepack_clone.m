function [nc]=ridgepack_clone(ncfile,newvar,newrec,newrecp,newendp,nearestdate)

% ridgepack_clone - Clones data and its metadata from a netCDF file to a MATLAB structure
%
% function [nc]=ridgepack_clone(ncfile,newvar,newrec,newrecp,newendp,nearestdate)
%
% INPUT:
%
% ncfile - Name of the netCDF file to be read. NOTE THAT ALL NETCDF 
%          INPUT FILES SHOULD END WITH THE SUFFIX ".nc".
%
% newvar - Name of a particular variable to be read.  This should
%          be included if you do not want to read in all data from the
%          netCDF file. It may be entered as either a string
%          specifying the single variable you wish to import, or
%          a cell array if you wish to import multiple variables.
%          Note that there is no need to specify 'supporting variables'
%          such as constant dimensional or coordinate data that help to describe
%          the variable you are importing if there is an unlimited dimension,
%          such as time, specified in the dataset.  These variables are 
%          automatically imported along with the variable you choose
%          so that it is fully self-described in the nc structure.
%          If you need to specify a particular hyperslab using newrec, but
%          wish to import all variables for that hyperslab that appear
%          in the netCDF file, enter newvar as an empty cell array: {}
%          (newvar is an optional argument)
%
% newrec - This input can have three possible types, each carrying
%          a particular meaning. If omitted, all of newvar is imported:
%          1) newrec is a number:
%          This specifies the particular time record to be read.
%          For example, newrec=3 will extract the third time record.
%          In this configuration, any time_bounds descriptors are NOT
%          processed.
%          2) newrec is a character string:
%          This specifies a particular time to be read. For example,
%          setting newrec='1-1-1997 00:00:00' will extract this date, 
%          or the nearest record in the netCDF file on the same day. 
%          This functionality only works where time is both a dimension
%          and variable in the netCDF file.  It does not work, for
%          example, for WRF model output files where time is specified
%          as a character string variable.
%          3) newrec is a cell array:
%          This specifies dimensions over which the extracted data 
%          are to be limited, and requires an additional argument
%          to specify the index of each slice to be provided in newrecp.
%          For example, to read a slice of data for a variable
%          with dimensions {'time','x','y','z'} on the y and z
%          plane, you would need to specifiy the values of time and
%          x for which the plane is defined, so newrec={'time','x'}.
%          The exact indices of time and x that define the plane are 
%          provided in newrecp For example newrecp={3,10} would specify 
%          for the 3rd and 10th indices for time and x, respectively. 
%          newrecp={3,'30'} would specify the 3rd index and nearest x  
%          to 30, and newrecp={'1-1-1997 00:00:0.0',10} would give the
%          the specified time and 10th x index. Note that when newrec
%          is left empty, newvar is extracted in complete form.
%          Note that for the latter to examples where actual values
%          are provided, the dimension must also be a variable in 
%          the netCDF file.
%          (newrec is an optional argument)
%
% Subsequent input for the special case of newrec being a cell array:
%
% newrecp - If newrec is a cell array, then newrecp is used to determine
%          the index for each of the dimensions in newrec for which a 
%          slice is being read from an ncfile for newvar, as described
%          for newrec.
%          (newrecp is an optional argument that should be specfied with newrec)
% 
% newendp - If newrec is a cell array, you may optionally choose to provide
%          an end bound on the netCDF data imported.  If newendp is omitted,
%          it is assumed that a single slice is being taken through a 
%          dimension.  However if you wish to extract a block of data, then
%	   newendp should be included to provide the end indices.  Hence
%          newrecp to newendp provides the slice required.  Note that 
%          newendp is NOT the 'count' often used in netCDF terms.  It
%          is the actual end index for which the data should be extracted.
%          (newrecp is an optional argument)
%
% nearestdate - This chooses the time unit ("year", "month" or "day") for 
%          which Icepack chooses the nearest possible time with the option 3
%          use of newrec. If you choose "year", then ridgepack_clone will take the 
%          nearest date in the same year.  If you choose "month", it will 
%          take data from the nearest date in the same month as specified, 
%          and if you choose "day", then the nearest 24 hours will be chosen. 
%          The default is "day".  If nearestdate is chosen as a number, then
%          this sets the nearest number of days to be chosen. The default is 
%          to only accept the nearest time within 24 hours of the specified
%          newrecp and newendp dates.
%
% If newvar and newrec are excluded, the entire netCDF file will be 
% extracted.  If newrec is excluded, only newvar will be extracted plus 
% any variables that describe its dimensions. If newrec is a cell
% array and newrecp is excluded, an error will result.
%
%
% OUTPUT:
%
% nc - netCDF data structure briefly described here:
%
% This function extracts data from a netCDF file 'ncfile' and 
% places it in a data structure "nc".  The data are 
% structured as they are within netCDF for variable 'name':
%
% nc.attributes         - global attributes
% nc.name.data          - data values
% nc.name.units         - units of the data
% nc.name.long_name     - long name of the data
% nc.name.fillvalue     - fill values substituted into data
%
% Where nc is the matlab data structure. Data with scale 
% factors and offsets are already adjusted upon being read 
% from the netCDF file. No further adjustment is required.  
% The missing_values and _fillvalues in the data are substituted 
% with NaNs (hence both are treated the same within Matlab).
% Further field names are possible and available, but are not
% shown here.  For more information on the netCDF structure, 
% see 'help ridgepack_struct'.
%
% Ridgepack Version 1.0
% Andrew Roberts, Naval Postgraduate School, March 2018 (afrobert@nps.edu)
%

% Set debug flag;
global debug;
if debug; disp(['Entering ',mfilename,'...']); end


% Run a check on input file name
if nargin>0 & ~ischar(ncfile)
 error('Netcdf file name is not a character variable');
elseif nnz(isstrprop(ncfile, 'alphanum'))==0;
 error('No name given for the netCDF file');
end


% Add .nc to filename if this has not already been done
ui=length(ncfile); li=max(1,ui-2);
if ui < 4 || not(strcmp(ncfile(li:ui),'.nc')) ; ncfile=[ncfile,'.nc']; end
disp(['Opening ',ncfile,' using netCDF V',netcdf.inqLibVers]);

% set default unlimited dimension name
oldtimename='time';

% ....................................................................
% Extract data from netCDF file 

% Open the netCDF file

try
 ncid=netcdf.open(ncfile,'NC_NOWRITE');
catch
 ncfileold=ncfile;
 ncfile=ncfile(1:end-3);
 try
  ncid=netcdf.open(ncfile,'NC_NOWRITE');
 catch
  xpwd=pwd;
  error(['Can''t open ',ncfileold,' in ',xpwd])
 end
end
 

% Get overall information about the netCDF file
[ndimss,nvars,natts,unlimdimid] = netcdf.inq(ncid);

% Get global attributes, assigning an appropriate attribute name
gid=netcdf.getConstant('NC_GLOBAL');
for i=1:natts;
 nattname = netcdf.inqAttName(ncid,gid,i-1);
 attname=nattname;

 % remove hyphen from attribute name if necessary
 dash=strfind(attname,'-');
 if ~isempty(dash)
  for i=1:length(dash); attname(dash(i):dash(i))='_'; end
 end

 % remove dot from attribute name if necessary
 dot=strfind(attname,'.');
 if ~isempty(dot)
  for i=1:length(dot); attname(dot(i):dot(i))='_'; end
 end

 % reduce the length of attname
 if length(attname)>60
  if debug; disp(['Reducing the length of the attribute ',attname]); end
  attname=attname(max(1,end-60):end);
 end

 % remove '_' from the beginning and end of each attribute
 underscore=strfind(attname,'_');
 if ~isempty(underscore) && underscore(1)==1; attname=attname(2:end); end
 if ~isempty(underscore) && any(underscore==length(attname)); attname=attname(1:end-1); end

 % check for replication
 if exist('nc') && isfield(nc.attributes,attname); attname=[attname,'_i']; end

 % assign the global attribute to the nc structure
 nc.attributes.(attname) = netcdf.getAtt(ncid,gid,nattname);
end

if nargin==1; % extract all data

 if debug; disp('Extracting all netCDF data from this file'); end

 newvar=cell(1,nvars);
 for i=1:nvars
  newvar{i}=netcdf.inqVar(ncid,i-1);
 end

else ; % get an individual record

 if debug; disp('Extracting only some selected data from this file'); end

 % set the variables to be extracted into a cell array if needed
 if ischar(newvar)
  newv{1}=newvar;
  newvar=newv;
 elseif ~iscell(newvar) 
  error('var must either be a character variable or cell array')
 elseif isempty(newvar)
  newvar=cell(1,nvars);
  for i=1:nvars
   newvar{i}=netcdf.inqVar(ncid,i-1);
  end
 end

end

% check input
if nargin==3 & ~(isnumeric(newrec) | ischar(newrec))
  error('newrec must be a numeric index or character variable in this configuration');
elseif nargin>=4 & ~iscell(newrec) 
  error('newrec must be a cell array in this configuration');
elseif nargin>=4 & ~iscell(newrecp) 
  error('newrecp must be a cell array in this configuration');
elseif nargin>=5 & ~iscell(newendp) 
  error('newendp must be a cell array in this configuration');
elseif nargin>6
  error('Too many arguments have been entered');
elseif nargin==3 & ~iscell(newrec) & isnumeric(newrec)
  quickpick=true;
else
  quickpick=false;
end

% append available added time descriptors to variables aside from time itself
% this is only done if the full dimension descriptors are used
if ~quickpick
 addedtimestuff={'time_bounds'};
 for k=1:length(addedtimestuff)
  for i=1:nvars
   if strcmpi(netcdf.inqVar(ncid,i-1),char(addedtimestuff{k}))
    newvar{length(newvar)+1}=char(addedtimestuff{k});
   end
  end
 end
end

% set default nearest date value to "day" 
daytolerance=1;
if nargin<6
 nearestdate='day'; 
elseif isnumeric(nearestdate)
 daytolerance=nearestdate;
 nearestdate='day'; 
elseif ~(strcmp(nearestdate,'year') | ...
         strcmp(nearestdate,'month') | ...
         strcmp(nearestdate,'day')) 
 error('nearestdate must be either ''year'',''month'', or ''day''')
end

% assign netCDF constant values to their names in a cell array
constnames={'NC_BYTE','NC_CHAR','NC_SHORT','NC_INT','NC_FLOAT','NC_DOUBLE'};
constvalues=zeros([length(constnames) 1]);
for cnam=1:length(constnames)
 constvalues(cnam)=netcdf.getConstant(char(constnames{cnam}));
end
[constvalues,I]=sort(constvalues);
constnames=constnames(I);

% set the count on the time variable
timecount=0;
timedims=[];
timevar=[];

% set to empty the list of dimension and variable IDs already loaded 
alldimids=[];
allvardid=[];

% This segment taken out from inside the mm loop, may have significance
% not understood when this was done.
if nargin>=3; rec=newrec; end
if nargin>=4; recp=newrecp; end
if nargin==4; endp=recp; end
if nargin>=5; endp=newendp; end

for mm=1:length(newvar)

  var=char(newvar{mm});

  % check that the variable exists
  try
   varid=netcdf.inqVarID(ncid,var);
  catch
   error(['Check variable name: Unable to find ''',var,''' in ',ncfile]);
  end

  % get dimension ids for the variable and dimension the nc.(var).dimension variable
  [varname,vtype,dimids,natts] = netcdf.inqVar(ncid,varid);
  dimids=dimids(end:-1:1); % reverse dimids so that the order is the same as for var data

  % check that all rec dimension variables apply to var
  if nargin>=4 && iscell(newrec); 
   for kkk=1:length(newrec)
    try
     did=netcdf.inqDimID(ncid,char(rec{kkk}));
    catch
     error(['ERROR: ',char(rec{kkk}),' Dimension specification error.']);
    end
   end
  end
 
  % compile list of all previously encountered dimensions
  dimids=setdiff(dimids,alldimids);
  alldimids=[alldimids setdiff(dimids,alldimids)];

  % Find dimension variables of var and extract them.  If a variable does not match the
  % name of the dimension, then an integer variable is added as a substitute.  This
  % potentially causes a problem if, for example, the dimension is called 'lon' and 
  % the associated variable is called 'longitude'. However this is likely to be remedied
  % by also importing all non-time dimensioned variables in the section following this.

  for m = 1:length(dimids);

    % find the variable that only has one dimid identical to one of dimids of var
    name=[];
    for jj=1:nvars
     [varname,xtype,did,ndatts] = netcdf.inqVar(ncid,jj-1);
     if length(did)==1 & did==dimids(m) 
      [dimname,dlength] = netcdf.inqDim(ncid,did);
      if strcmpi(dimname,varname); name=varname; break; end
     end
    end

    if isempty(name) 

     [name,dlength] = netcdf.inqDim(ncid,dimids(m));

     % check name for leading non-alphabetic characters
     charlog=isstrprop(name,'alpha');
     for i=1:length(name)
      if charlog(i); modifiedname=name(i:end); break; end
      if i==length(name); modifiedname=['dim',num2str(dimids(m))]; end
     end
     name=modifiedname;
     charlog=isstrprop(name,'alphanum');
     name(not(charlog))='_';

     if debug; disp(['Creating ',name,' as a dimension variable']); end;

     % add incidentals to the nc name
     nc.(name).long_name=name;
     nc.(name).type='NC_INT';
     nc.(name).dimension={name};

     % construct dimensional array 
     if nargin>=4

      recindex=find(strcmp(name,rec));
      if length(recindex)>1
       error(['More than one match for ',name])
      elseif isempty(recindex) 
       nc.(name).data=[1:dlength];
       if debug; disp(['Setting indexed dimension ',...
                       name,' from 1 to ',num2str(dlength)]); end
      elseif ~isnumeric(recp{recindex}) | ~isnumeric(endp{recindex})
       if (ischar(recp{recindex}) | ischar(endp{recindex}))
        disp(['Date search not applicable to ',ncfile,'.'])
        if nargin==4
         error(['Change ''',char(recp{recindex}),''' to a numeric index.'])
        elseif nargin>4
         error(['Change ''',char(recp{recindex}),''' & ',...
                       char(endp{recindex}),' to numeric indices.'])
        end
       else
        error('Processing recp and endp.')
       end
      else 
       if recp{recindex}>0 & endp{recindex}<=dlength & recp{recindex}<=endp{recindex}
        nc.(name).data=[recp{recindex}:1:endp{recindex}];
        if debug; disp(['Setting indexed dimension ',name,...
                       ' from ',num2str(recp{recindex}),...
                       ' to ',num2str(endp{recindex})]); end
       else
        error(['Dimension ',name,' length is ',num2str(dlength)])
       end
      end

     % if a time record number is specified for rec, assign it to recp
     elseif quickpick & strcmpi(name,'time') 

      if isnumeric(newrec)
       recp{m}=newrec;
       endp{m}=newrec;
       if ~iscell(rec); clear rec; end
       rec{m}=name;
       nc.(name).data=recp{m};
       if debug; disp(['Selecting indexed time dimension ',num2str(recp{1})]); end
      elseif ischar(newrec)
       disp(['Date search not applicable to ',ncfile,'.']) 
       error(['Change ''',rec,''' to a numeric index.'])
      else
       error('rec specification error.')
      end

     % you landed here if the third argument specifies a record number, but a
     % requested variable (var) has a variable-less dimension other then time
     elseif quickpick

      if isnumeric(newrec)
       recp{m}=1;
       endp{m}=dlength;
       if ~iscell(rec); clear rec; end
       rec{m}=name;
       nc.(name).data=[1:dlength]';
       if debug; 
         disp(['Selecting dimension ',name,' from with length ',num2str(endp{1})]); 
       end
      else
       error('newrec specification error.')
      end

     % otherwise assign full length vector for dimension
     else

      nc.(name).data=[1:dlength];
      if debug; 
       disp(['Setting indexed dimension ',name,' from 1 to ',num2str(dlength)]); 
      end

     end

     % rename time variable to 'time'
     if strcmpi(name,'time') & ~isfield(nc,oldtimename); 
      if(debug);disp('Time dimension not name-matched with a variable');end
      oldtimename=name;
     end

    else

     if xtype==netcdf.getConstant('NC_CHAR');
      error(['Unable to assign a dimension to the character dimension ',name])
     elseif debug
      disp(['Extracting ',name,' as a dimension variable']);
     end

     % get the variable ID of name
     vardid=netcdf.inqVarID(ncid,name);

     % check name for leading non-alphabetic characters
     charlog=isstrprop(name,'alpha');
     for i=1:length(name)
      if charlog(i); modifiedname=name(i:end); break; end
      if i==length(name); modifiedname=['var',num2str(vardid)]; end
     end
     name=modifiedname;
     charlog=isstrprop(name,'alphanum');
     name(not(charlog))='_';
     if debug; disp(['Modified name is ',name]); end;

     % compile list of all variables that have already been loaded
     allvardid=[allvardid setdiff(vardid,allvardid)];

     % get attributes of name
     for attnum=1:ndatts
      attname=netcdf.inqAttName(ncid,vardid,attnum-1);
      % check attname for leading non-alphabetic characters
      charlog=isstrprop(attname,'alpha');
      for i=1:length(attname)
       if charlog(i); modifiedattname=attname(i:end); break; end
       if i==length(attname); modifiedattname=['att',num2str(attnum)]; end
      end
      charlog=isstrprop(modifiedattname,'alphanum');
      modifiedattname(not(charlog))='_';
      nc.(name).(modifiedattname)=netcdf.getAtt(ncid,vardid,attname);
     end
     nc.(name).type=constnames{xtype};
     nc.(name).dimension={name};

     % get the time coordinate
     if unlimdimid==did & ...
        isfield(nc.(name),'units') & ...
        ~isempty(strfind(name,'time')) & ...
        (~isempty(strfind(nc.(name).units,'hours since')) | ...
         ~isempty(strfind(nc.(name).units,'Hours since')) | ...
         ~isempty(strfind(nc.(name).units,'HOURS SINCE')) | ...
         ~isempty(strfind(nc.(name).units,'minutes since')) | ...
         ~isempty(strfind(nc.(name).units,'Minutes since')) | ...
         ~isempty(strfind(nc.(name).units,'MINUTES SINCE')) | ...
         ~isempty(strfind(nc.(name).units,'seconds since')) | ...
         ~isempty(strfind(nc.(name).units,'Seconds since')) | ...
         ~isempty(strfind(nc.(name).units,'SECONDS SINCE')) | ...
         ~isempty(strfind(nc.(name).units,'days since')) | ...
         ~isempty(strfind(nc.(name).units,'Days since')) | ...
         ~isempty(strfind(nc.(name).units,'DAYS SINCE')));

      if(debug); 
       disp(['Working on unlimited time variable called ',name]); 
      end

      % register the original name of the time dimension variable
      oldtimename=name;

      % set time dimension and variable id
      timecount=timecount+1;
      timevar(timecount)=vardid;
      timedims(timecount)=did;
 
      % If rec is not a structure, convert a character rec variable 
      % to an actual time and index
      if nargin==3 && ischar(newrec) ; 

       index=newrec;
       if ~iscell(rec); clear rec; end
       rec{m}=name;
       recp{m}=index;
       endp{m}=index;
 
      % OR convert a specific record request to the rec structure
      elseif quickpick ;

       if newrec>dlength | newrec<1;
        error([dimname,' dimension has limits 1 to ',num2str(dlength),' in ',ncfile])
       end
       index=newrec;
       if ~iscell(rec); clear rec; end
       rec{m}=name;
       recp{m}=index;
       endp{m}=index;

      % OR just set endp if it has not been assigned
      elseif nargin==4

       endp=recp;

      end

      % Convert a time string in the structure to an index
      if nargin>2 && iscell(rec) ;

       if(debug); disp(['Converting time string to an index for ',name]); end

       index=find(strcmp(rec,name));

       % set value of recp
       if ~isempty(index) && ischar(recp{index})

        numberrecs=num2str(dlength);

        if debug; 
          disp(['Searching ',numberrecs,' for matching ',char(recp{index}),'.']); 
        end

        try
         timematch=datenum(recp{index});
        catch
         error(['Unable to process the date string ''',char(recp{index}),'''']);
        end

        input_data=netcdf.getVar(ncid,vardid);
        try; 
         calendar=netcdf.getAtt(ncid,vardid,'calendar'); 
        catch;
         calendar='standard'; 
        end
        try;
         units=netcdf.getAtt(ncid,vardid,'units'); 
        catch; 
         error(['No units for ',name]); 
        end
        try;
         scale_factor=netcdf.getAtt(ncid,vardid,'scale_factor');
        catch; 
         scale_factor=1; 
        end
        try; 
         add_offset=netcdf.getAtt(ncid,vardid,'add_offset'); 
        catch;
         add_offset=0;
        end
        try; 
         fillvalue=netcdf.getAtt(ncid,vardid,'_fillvalue'); 
        catch;
         fillvalue=NaN; 
        end
        if isnan(fillvalue)
         try; 
          fillvalue=netcdf.getAtt(ncid,vardid,'_FillValue'); 
         catch; 
          fillvalue=NaN; 
         end
        end
        try; 
          missing_value=netcdf.getAtt(ncid,vardid,'missing_value'); 
        catch; 
          missing_value=NaN; 
        end
        try
         timeshape=size(double(input_data(:))*scale_factor+add_offset);
         timevec=ridgepack_timeconvert(double(input_data(:))*scale_factor+add_offset,...
                               units,calendar,0);
         timevec=reshape(timevec,timeshape);
         timevec(input_data==fillvalue | input_data==missing_value )=NaN;
        catch
	 error(['Unable to convert ',name]);
        end

        indexp=find(abs(timevec-timematch)==min(abs(timevec-timematch)));
        indexp=indexp(1); % in case time match is half way between two values
        dvecnearest=datevec(timevec(indexp));
        dvecmatch=datevec(timematch);
        if (strcmp(nearestdate,'year') & ...
            dvecnearest(1)==dvecmatch(1)) | ...
           (strcmp(nearestdate,'month') & ...
            dvecnearest(1)==dvecmatch(1) & ...
            dvecnearest(2)==dvecmatch(2)) | ...
           (strcmp(nearestdate,'day') & ...
            abs(datenum(dvecnearest)-datenum(dvecmatch))<=daytolerance)
         disp(['Nearest to ',char(recp{index}),' is ',datestr(timevec(indexp))]);
        else
         warning('Unable to find the date requested')
         error(['Nearest to ',char(recp{index}),' is ',datestr(timevec(indexp))]);
        end
        
        oldrecp=char(recp{index});
        recp{index}=indexp;

        if nargin<5 | strcmp(oldrecp,char(endp{index}))
         endp{index}=recp{index};
        elseif ~ischar(endp{index});
	 error('newendp must have the same class as newrecp for each dimension')
	else
         numberrecs=num2str(dlength);

         if debug;disp(['Searching ',numberrecs,' records for matching endp time.']); end

         timematch=datenum(endp{index});
         indexp=find(abs(timevec-timematch)==min(abs(timevec-timematch)));
         indexp=indexp(end); % in case time match is half way between two values
         dvecnearest=datevec(timevec(indexp));
         dvecmatch=datevec(timematch);
         if (strcmp(nearestdate,'year') & ...
             dvecnearest(1)==dvecmatch(1)) | ...
            (strcmp(nearestdate,'month') & ...
             dvecnearest(1)==dvecmatch(1) & ...
             dvecnearest(2)==dvecmatch(2)) | ...
            (strcmp(nearestdate,'day') & ...
             abs(datenum(dvecnearest)-datenum(dvecmatch))<=daytolerance)
          disp(['Nearest to ',char(endp{index}),' is ',datestr(timevec(indexp))]);
         else
          warning('Unable to find the date requested')
          error(['Nearest to ',char(endp{index}),' is ',datestr(timevec(indexp))]);
         end

         endp{index}=indexp;

        end

       end
 
       if ~isempty(index) && (length(recp{index})>1 | length(endp{index})>1)
	error('Specify only one integer of the start and end index of each dimension')
       elseif ~isempty(index) & recp{index}>endp{index}
	tempp=recp{index};
	recp{index}=endp{index};
	endp{index}=tempp;
       end

      elseif nargin>2
       error('programming error: rec not converted to a structure');
      end

      % get entire time vectors if rec is not specified or time is not in rec
      if nargin<=2 | (nargin>2 & isempty(index))

        if debug; disp('Getting entire time vector'); end

        try
         input_data=netcdf.getVar(ncid,vardid);
         try;
          scale_factor=netcdf.getAtt(ncid,vardid,'scale_factor');
         catch;
           scale_factor=1;
         end
         try;
          add_offset=netcdf.getAtt(ncid,vardid,'add_offset');
         catch;
         add_offset=0;
         end
         try;
          fillvalue=netcdf.getAtt(ncid,vardid,'_fillvalue');
         catch;
          fillvalue=NaN;
         end
         if isnan(fillvalue)
          try;
           fillvalue=netcdf.getAtt(ncid,vardid,'_FillValue');
          catch;
           fillvalue=NaN;
          end
	 end
         try;
          missing_value=netcdf.getAtt(ncid,vardid,'missing_value'); 
         catch;
          missing_value=NaN; 
         end
         nc.(name).data=double(input_data)*scale_factor+add_offset;
         nc.(name).data(input_data==fillvalue | input_data==missing_value )=NaN;
        catch
         error(['Unable to extract the indexed data for ',name])
        end
 
      % OR just get the indexed time
      elseif nargin>2 & ~isempty(index)

        if debug; disp('Getting indexed time'); end

        timeischar=false;
        try
         input_data=netcdf.getVar(ncid,vardid,recp{index}-1,endp{index}-recp{index}+1);
         try;
          scale_factor=netcdf.getAtt(ncid,vardid,'scale_factor'); 
         catch;
          scale_factor=1; 
         end
         try;
          add_offset=netcdf.getAtt(ncid,vardid,'add_offset'); 
         catch;
          add_offset=0; 
         end
         try;
          fillvalue=netcdf.getAtt(ncid,vardid,'_fillvalue'); 
         catch;
          fillvalue=NaN; 
         end
         if isnan(fillvalue)
          try;
           fillvalue=netcdf.getAtt(ncid,vardid,'_FillValue'); 
          catch;
           fillvalue=NaN; 
          end
         end
         try;
          missing_value=netcdf.getAtt(ncid,vardid,'missing_value'); 
         catch;
          missing_value=NaN; 
         end
         nc.(name).data=double(input_data)*scale_factor+add_offset;
         nc.(name).data(input_data==fillvalue | input_data==missing_value )=NaN;
        catch
         error(['Unable to extract the indexed data for ',name])
        end
	nc.(name).data=double(nc.(name).data);

      % error message if reached here without satisfaction
      else

       error('Inputs are incorrect')

      end 

      % convert time variable to small case time, and keep previous name of variable
      if strcmpi(name,'time')&~strcmp(name,'time'); 
	nc.time=nc.(name); 
	nc=rmfield(nc,name); 
      end
 
     % for non-time dimensions, follow a similar procedure
     else

      if(debug);
       disp(['Working on dimensions not recognised as time: ',name]); 
      end

      if nargin==4
       endp=recp;
      end

      % dimension variables for selected grid points
      if nargin>2 && iscell(rec)


       index=find(strcmp(rec,name));

       if ~isempty(index) && ischar(recp{index})

        pointmatch=str2num(recp{index});

	nc.(name).data=netcdf.getVar(ncid,vardid);
        try;
         scale_factor=netcdf.getAtt(ncid,vardid,'scale_factor'); 
        catch;
         scale_factor=1; 
        end
        try;
         add_offset=netcdf.getAtt(ncid,vardid,'add_offset'); 
        catch;
         add_offset=0; 
        end
 	pointvec = double(nc.(name).data) * scale_factor + add_offset;

        indexp=find(abs(pointvec-pointmatch)==min(abs(pointvec-pointmatch)));
        disp(['Nearest ',rec{index},' to ',recp{index},' is ',num2str(pointvec(indexp))]);
        recp{index}=indexp;

        if nargin<5

         endp{index}=recp{index};

        elseif ~ischar(endp{index});

	 error('newendp must have the same class as newrecp for each dimension')

	else

         pointmatch=str2num(endp{index});

	 nc.(name).data=netcdf.getVar(ncid,vardid);
         try;
          scale_factor=netcdf.getAtt(ncid,vardid,'scale_factor'); 
         catch;
          scale_factor=1; 
         end
         try;
          add_offset=netcdf.getAtt(ncid,vardid,'add_offset'); 
         catch;
          add_offset=0; 
         end
 	 pointvec = double(nc.(name).data) * scale_factor + add_offset;

         indexp=find(abs(pointvec-pointmatch)==min(abs(pointvec-pointmatch)));
         disp(['Nearest ',rec{index},' to ',endp{index},' is ',num2str(pointvec(indexp))]);
         endp{index}=indexp;

        end

       end

       if ~isempty(index) && recp{index}>endp{index}
	 tempp=recp{index};
	 recp{index}=endp{index};
	 endp{index}=tempp;
       end


      elseif nargin>2;

       index=[];

      end

      % get entire vector if rec is not specified or name is not in rec
      if nargin<=2 | (nargin>2 & isempty(index))

        if(debug); 
          disp(['Getting entire ',name,' variable (rec not specified or name not in rec)']);
        end

        try
	 nc.(name).data=netcdf.getVar(ncid,vardid);
        catch
         error(['Can''t find ',name]);
        end
        try;
         scale_factor=netcdf.getAtt(ncid,vardid,'scale_factor'); 
        catch;
         scale_factor=1; 
        end
        try;
         add_offset=netcdf.getAtt(ncid,vardid,'add_offset'); 
        catch;
         add_offset=0; 
        end
 	nc.(name).data = double(nc.(name).data) * scale_factor + add_offset;

      % OR just get the chosen value
      elseif nargin>2 && ~isempty(index)

        if(debug);
         disp(['Getting entire ',name,' variable (rec not specified or name not in rec)']);
        end

        try
	 nc.(name).data=netcdf.getVar(ncid,vardid,recp{index}-1,endp{index}-recp{index}+1);
        catch
         error(['Can''t find ',name,': ',num2str(recp{index}-1),...
                ' to ',num2str(endp{index})]);
        end
        try;
         scale_factor=netcdf.getAtt(ncid,vardid,'scale_factor');
        catch;
         scale_factor=1; 
        end
        try;
         add_offset=netcdf.getAtt(ncid,vardid,'add_offset'); 
        catch;
         add_offset=0; 
        end
 	nc.(name).data = double(nc.(name).data) * scale_factor + add_offset;
 
      % error message if reached here without satisfaction
      else

       error('Inputs are incorrect')

      end 
    
     end

    end

    moddimname{dimids(m)+1}=name;

  end

  % Extract requested variable 'var' and variables that do not have an unlimited dimension 
  % and share the same other dimensions with 'var'. This is done so that useful variables 
  % related to var, such as masks or coordinate variables are sucked into the nc structure, 
  % even if they are not explicitly mentioned in the attributes of var.  This only occurs
  % if 'var' has a time dimension.

  for i = 1:nvars;

   vardid=i-1;
 
   [name,vartype,dimids,nvatts]=netcdf.inqVar(ncid,vardid);

   if  ~isempty(setdiff(vardid,allvardid)) & ...
        length(dimids)==length(intersect(dimids,alldimids)) & ...
        (strcmp(name,var) | strcmp(name,'time_bounds') | ...
        (isempty(intersect(dimids,[timedims unlimdimid])) & ...
        unlimdimid>-1 & ~isempty(timedims))) ;

    % add to list of variables already encountered
    allvardid=[allvardid setdiff(vardid,allvardid)];

    % check name for leading non-alphabetic characters
    charlog=isstrprop(name,'alpha');
    for j=1:length(name)
     if charlog(j); modifiedname=name(j:end); break; end
     if j==length(name); modifiedname=['var',num2str(vardid)]; end
    end
    name=modifiedname;
    charlog=isstrprop(name,'alphanum');
    name(not(charlog))='_';

    % get attributes of name, removing leading non-alphabetic characters from attname
    for attnum=1:nvatts
     attname=netcdf.inqAttName(ncid,vardid,attnum-1);
     charlog=isstrprop(attname,'alpha');
     for j=1:length(attname)
      if charlog(j); modifiedattname=attname(j:end); break; end
      if j==length(attname); modifiedattname=['att',num2str(attnum)]; end
     end
     nc.(name).(modifiedattname)=netcdf.getAtt(ncid,vardid,attname);
    end
    nc.(name).type=constnames{vartype};

    % assign dimension cell array
    for dk=1:length(dimids)
     nc.(name).dimension{dk}=char(moddimname{dimids(dk)+1});
    end

    % check that listed coordinate variables do not have time as a dimension,
    % but if they do, get them at the same time as getting the variable. Also,
    % change the name of the coordinate listing to be consistent with changes
    % in variable names (removal of nonalphabetic characters from name).
    getcoords=0;

    if isfield(nc.(name),'coordinates')
     cords=char(nc.(name).coordinates);
     if ~isempty(cords)

      % anticpate getting coordinate variables
      getcoords=1;

      % extract individual coordinate variables from the coordinate list
      for chpos=1:10
       [coords{chpos},cords]=strtok(cords,' ');
       % Th line below was changed because lats and longs were going missing 
       % and ** was added below
       %if strcmp(cords,''); coords=coords(1:chpos-1); break; end
       if strcmp(cords,''); coords=coords(1:chpos); break; end
      end

      % ** remove blanks from any part of coords
      cords=coords;
      k=0;
      for chpos=1:length(cords)
       if ~strcmp(char(cords{chpos}),' ')
        k=k+1;
        coords{k}=cords{chpos};
       end
      end

      % now interrogate netCDF file for coordinate variables
      for itc=1:length(coords)

       % extract coordinate information
       try
        coordid(itc)=netcdf.inqVarID(ncid,char(coords{itc}));
        [cname,cvartype,cdimids,cnvatts]=netcdf.inqVar(ncid,coordid(itc));
       catch
	disp(['Cannot locate coordinate ',char(coords{itc}),' for variable ',name])
       end

       % Check for time coords in dimensions.  If it exists, get coordinate variables
       % at the same time as getting the variable because it won't be picked up unless
       % it is explicitly listed.  Also, check that the dimensions for the listed
       % coordinate variables are the same as for the variable.  If not, do not
       % implicitly pick up the coordinate variable.
       if ~isempty(intersect(cdimids,[timedims unlimdimid]))
	for ic=1:length(cdimids)
         dimfind=find(dimids==cdimids(ic));
	 if isempty(dimfind);
	  disp(['Unable to find a dimension for ',cname,' in ',name,...
                ' dimensions, therefore'])
	  disp(['coordinate variables for ',name,...
                ' must be explicitly listed in ridgepack_clone list.'])
          getcoords=0;
	 end
	end 
       else
        getcoords=0;
       end

       % While in the loop, rename the listed coordinate to that with a corrected
       % name for the listed variable.  This needs to be done regardless of whether
       % or not the coordinate is time dependent
       modifiedname=cname;
       charlog=isstrprop(cname,'alpha');
       for j=1:length(cname)
        if charlog(j); modifiedname=cname(j:end); break; end
        if j==length(cname); modifiedname=['var',num2str(coordid(itc))]; end
       end
       newcoords{itc}=modifiedname;

      end


      % rewrite the corrected coordinate strings
      if ~exist('newcoords','var') | isempty(newcoords) 
       nc.(name)=rmfield(nc.(name),'coordinates');
       getcoords=false;
       disp(['Coordinates attribute broken for ',name])
      else
       nc.(name).coordinates=ridgepack_cellcat(newcoords);
       nc.(name).coordinates=nc.(name).coordinates(2:end);
      end


     else

      disp('(The listed coordinates attribute appears to be empty for ',name,')')

     end

     % If picking up coordinate variables affiliated with the variable itself, 
     % add the coordinates to the list of variables already encountered.
     % Then set up coordinate variable

     if getcoords; 

       allvardid=[allvardid setdiff(coordid,allvardid)];

       % get dimensions for time-dependent coordinate variable
       for itc=1:length(coordid)

        cvardid=coordid(itc);

        [cname,cvartype,cdimids,cnvatts]=netcdf.inqVar(ncid,cvardid);
	ccname{itc}=cname;

        if debug; disp(['Setting up ',cname,' coordinate variable']); end;

        % check cname for leading non-alphabetic characters
        charlog=isstrprop(cname,'alpha');
        for j=1:length(cname)
         if charlog(j); modifiedname=cname(j:end); break; end
         if j==length(cname); modifiedname=['var',num2str(cvardid)]; end
        end
        cname=modifiedname;
        charlog=isstrprop(cname,'alphanum');
        cname(not(charlog))='_';

        % get attributes of cname, removing leading non-alphabetic characters from cattname
        for attnum=1:cnvatts
         cattname=netcdf.inqAttName(ncid,cvardid,attnum-1);
         charlog=isstrprop(cattname,'alpha');
         for j=1:length(cattname)
          if charlog(j); cmodifiedattname=cattname(j:end); break; end
          if j==length(cattname); cmodifiedattname=['att',num2str(attnum)]; end
         end
         charlog=isstrprop(cmodifiedattname,'alphanum');
         cmodifiedattname(not(charlog))='_';
         nc.(cname).(cmodifiedattname)=netcdf.getAtt(ncid,cvardid,cattname);
        end

        % rassign coordinate string to the same as name where there have been changes
        nc.(cname).coordinates=nc.(name).coordinates;

        % assign coordinate variable type
	nc.(cname).type=constnames{cvartype};

        % assign coordinate dimension
        for dk=1:length(cdimids)
         nc.(cname).dimension{dk}=char(moddimname{cdimids(dk)+1});
        end

       end

     else

	clear coordst

     end

    end % end of time-dependent coodinate variable checking

    if  vartype~=netcdf.getConstant('NC_CHAR');

     % Get data selection of non-character variables
     if nargin>2 & ~isempty(dimids)

      start=zeros(size(dimids));
      count=zeros(size(dimids));
      for jj=1:length(dimids)
       [dimname,count(jj)] = netcdf.inqDim(ncid,dimids(jj));
       for kk=1:length(rec)
        if strcmpi(rec{kk},dimname)
         start(jj)=recp{kk}-1;
  	 count(jj)=endp{kk}-recp{kk}+1;
        end
       end
      end
 
      if debug;
        disp(['Extracting ',name,' as a variable: start=[',num2str(start),'] count=[',...
                      num2str(count),']']);
      end

      try
       input_data=netcdf.getVar(ncid,vardid,start,count);
       try; 
        scale_factor=netcdf.getAtt(ncid,vardid,'scale_factor'); 
       catch;
        scale_factor=1; 
       end
       try;
        add_offset=netcdf.getAtt(ncid,vardid,'add_offset'); 
       catch;
        add_offset=0; 
       end
       try;
        fillvalue=netcdf.getAtt(ncid,vardid,'_fillvalue'); 
       catch;
        fillvalue=NaN; 
       end
       if isnan(fillvalue)
        try;
         fillvalue=netcdf.getAtt(ncid,vardid,'_FillValue'); 
        catch;
         fillvalue=NaN; 
        end
       end
       try;
        missing_value=netcdf.getAtt(ncid,vardid,'missing_value'); 
       catch;
        missing_value=NaN; 
       end
       nc.(name).data = double(input_data) * scale_factor + add_offset;
       nc.(name).data(input_data==fillvalue | input_data==missing_value )=NaN;

      catch

       input_data=netcdf.getVar(ncid,vardid,start,count);
       error(['Unable to load ',name,', probably out of memory']);

      end

      if getcoords

       allvardid=[allvardid setdiff(coordid,allvardid)];

       % get dimensions for time-dependent coordinate variable
       for itc=1:length(coordid)

        cvardid=coordid(itc);

        [cname,cvartype,cdimids,cnvatts]=netcdf.inqVar(ncid,cvardid);
	ccname{itc}=cname;

        % check cname for leading non-alphabetic characters
        charlog=isstrprop(cname,'alpha');
        for j=1:length(cname)
         if charlog(j); modifiedname=cname(j:end); break; end
         if j==length(cname); modifiedname=['var',num2str(cvardid)]; end
        end
        cname=modifiedname;
        charlog=isstrprop(cname,'alphanum');
        cname(not(charlog))='_';

        % get attributes of cname, removing leading non-alphabetic characters from cattname
        for attnum=1:cnvatts
         cattname=netcdf.inqAttName(ncid,cvardid,attnum-1);
         charlog=isstrprop(cattname,'alpha');
         for j=1:length(cattname)
          if charlog(j); cmodifiedattname=cattname(j:end); break; end
          if j==length(cattname); cmodifiedattname=['att',num2str(attnum)]; end
         end
         charlog=isstrprop(cmodifiedattname,'alphanum');
         cmodifiedattname(not(charlog))='_';
         nc.(cname).(cmodifiedattname)=netcdf.getAtt(ncid,cvardid,cattname);
        end

        % rassign coordinate string to the same as name where there have been changes
        nc.(cname).coordinates=nc.(name).coordinates;

        % assign coordinate variable type
	nc.(cname).type=constnames{cvartype};

        % assign coordinate dimension
        for dk=1:length(cdimids)
         nc.(cname).dimension{dk}=char(moddimname{cdimids(dk)+1});
        end

        % determine bounds on coordinate data, only reading in one time slice
        cstart=zeros(size(cdimids));
        ccount=zeros(size(cdimids));
        for jj=1:length(cdimids)
         [cdimname,ccount(jj)] = netcdf.inqDim(ncid,cdimids(jj));
         if strcmpi(cdimname,'time')
          cstart(jj)=0;
  	  ccount(jj)=1;
         else
          for kk=1:length(rec)
           if strcmpi(rec{kk},cdimname)
            cstart(jj)=recp{kk}-1;
  	    ccount(jj)=endp{kk}-recp{kk}+1;
           end
          end
         end
        end

        if debug; 
          disp(['Extracting ',cname,' as a coordinate: cstart=[',num2str(cstart),...
                        '] ccount=[',num2str(ccount),']']);
        end

        try
         input_data=netcdf.getVar(ncid,cvardid,cstart,ccount);
         try; 
          scale_factor=netcdf.getAtt(ncid,cvardid,'scale_factor'); 
         catch;
           scale_factor=1; 
         end
         try;
           add_offset=netcdf.getAtt(ncid,cvardid,'add_offset'); 
         catch;
           add_offset=0; 
         end
         try;
          fillvalue=netcdf.getAtt(ncid,cvardid,'_fillvalue');
         catch; 
          fillvalue=NaN;
         end
         if isnan(fillvalue)
          try;
           fillvalue=netcdf.getAtt(ncid,cvardid,'_FillValue'); 
          catch;
           fillvalue=NaN; 
          end
         end
         try;
          missing_value=netcdf.getAtt(ncid,cvardid,'missing_value'); 
         catch;
          missing_value=NaN; 
         end
         nc.(cname).data = double(input_data) * scale_factor + add_offset;
         nc.(cname).data(input_data==fillvalue | input_data==missing_value )=NaN;
        catch
         error(['Unable to load coordinate variable ',cname]);
        end

        % Remove time variable from coordinate dimension if it exists
        tindex=find(strcmpi('time',nc.(cname).dimension));
        if ~isempty(tindex) & length(nc.(cname).dimension)>1
         xorder=[1:ndims(nc.(cname).data)];
         yorder=[xorder(xorder~=tindex) tindex];
         nc.(cname).dimension=nc.(cname).dimension(xorder(xorder~=tindex));
        end

       end

      end

     % OR get the entire dataset for var
     else
 
      try

       if debug; disp(['Extracting entire ',name,' record as a variable']); end;

       input_data=netcdf.getVar(ncid,vardid);
       try;
        scale_factor=netcdf.getAtt(ncid,vardid,'scale_factor'); 
       catch;
        scale_factor=1; 
       end
       try;
        add_offset=netcdf.getAtt(ncid,vardid,'add_offset'); 
       catch;
        add_offset=0; 
       end
       try;
        fillvalue=netcdf.getAtt(ncid,vardid,'_fillvalue'); 
       catch;
        fillvalue=NaN; 
       end
       if isnan(fillvalue)
        try;
         fillvalue=netcdf.getAtt(ncid,vardid,'_FillValue');
        catch;
         fillvalue=NaN;
        end
       end
       try;
        missing_value=netcdf.getAtt(ncid,vardid,'missing_value'); 
       catch;
        missing_value=NaN; 
       end
       nc.(name).data = double(input_data) * scale_factor + add_offset;
       nc.(name).data(input_data==fillvalue | input_data==missing_value )=NaN;

       if debug; 
        nc.(name)
       end

      catch

       error(['Unable to load ',name,', probably out of memory']);

      end

      if getcoords

       for itc=1:length(coordid)

        cname=char(ccname{itc});
        cvardid=coordid(itc);

        if debug; disp(['Extracting entire ',cname,' coordinate variable']); end;

        try

         input_data=netcdf.getVar(ncid,cvardid);
         try;
          scale_factor=netcdf.getAtt(ncid,cvardid,'scale_factor'); 
         catch; 
          scale_factor=1; 
         end
         try;
          add_offset=netcdf.getAtt(ncid,cvardid,'add_offset'); 
         catch; 
          add_offset=0; 
         end
         try;
          fillvalue=netcdf.getAtt(ncid,cvardid,'_fillvalue');
         catch; 
          fillvalue=NaN; 
         end
         if isnan(fillvalue)
          try;
           fillvalue=netcdf.getAtt(ncid,cvardid,'_FillValue');
          catch;
           fillvalue=NaN;
          end
         end
         try; 
          missing_value=netcdf.getAtt(ncid,cvardid,'missing_value'); 
         catch; 
          missing_value=NaN; 
         end
         nc.(cname).data = double(input_data) * scale_factor + add_offset;
         nc.(cname).data(input_data==fillvalue | input_data==missing_value )=NaN;

        catch

         error(['Unable to load coordinate variable ',cname]);

        end

       end
      
      end


     end

    elseif ~isempty(dimids) % retrieve character variables

     if debug; disp(['Extracting ',name,' as a character variable.']); end
     try
      nc.(name).data=netcdf.getVar(ncid,vardid);
     catch
      error(['Unable to load the character variable ',name]);
     end

    end

   end

  end

end


% ....................................................................
% Get variable names
[variablenames,numbervariables]=ridgepack_name(nc);
disp(['Importing ',num2str(numbervariables),' variable(s) into the netCDF structure']);

% ....................................................................
% Check for limited dimension variable called time and replace name with timexxxx
replacetime=0;
for i=1:numbervariables
 name=char(variablenames{i});
 if strcmp(name,'time') & ~strcmp(name,oldtimename)
  disp('Replacing time variable with timexxxx')
  nc.timexxxx=nc.time; nc=rmfield(nc,'time');
  replacetime=1;
 end
end

% replace limited dimension time with substitute name timexxxx for each variable
if replacetime
 [variablenames,numbervariables]=ridgepack_name(nc);
 for i=1:numbervariables
  name=char(variablenames{i});
  if any(strcmp(nc.(name).dimension,'time'))
   nc.(name).dimension{findstrcmp(nc.(name).dimension,'time')}='timexxxx';
  end
 end
end


% ....................................................................
% Turn time coordinates into serial matlab time, including special cases
for i=1:numbervariables

 name=char(variablenames{i});

 % read standard records
 if isfield(nc.(name),'units') & ...
    isfield(nc.(name),'data') & ...
    ~isempty(strfind(name,'time')) & ...
    (~isempty(strfind(nc.(name).units,'hours since')) | ...
     ~isempty(strfind(nc.(name).units,'Hours since')) | ...
     ~isempty(strfind(nc.(name).units,'HOURS SINCE')) | ...
     ~isempty(strfind(nc.(name).units,'minutes since')) | ...
     ~isempty(strfind(nc.(name).units,'Minutes since')) | ...
     ~isempty(strfind(nc.(name).units,'MINUTES SINCE')) | ...
     ~isempty(strfind(nc.(name).units,'seconds since')) | ...
     ~isempty(strfind(nc.(name).units,'Seconds since')) | ...
     ~isempty(strfind(nc.(name).units,'SECONDS SINCE')) | ...
     ~isempty(strfind(nc.(name).units,'days since')) | ...
     ~isempty(strfind(nc.(name).units,'Days since')) | ...
     ~isempty(strfind(nc.(name).units,'DAYS SINCE')));

  if (debug); disp('Assigning calendar to time dimension'); end

  if isempty(nc.(name).data); error([name,' is empty']); end

  if ~isfield(nc.(name),'calendar') & any(strcmpi('time',nc.(name).dimension)) & ...
     ~strcmp(name,'time') & isfield(nc,'time') & isfield(nc.time,'calendar'); 
   if debug; disp(['Setting ',name,' calendar to ',char(nc.time.calendar)]); end
   nc.(name).calendar=nc.time.calendar;
  elseif ~isfield(nc.(name),'calendar'); 
   if debug; disp(['Setting ',name,' calendar to gregorian']); end
   nc.(name).calendar='gregorian'; 
  end

  if not(isstrprop(nc.(name).units,'alphanum'));
    disp(char(nc.(name).units));
    error('units are not alphanumeric for time');
  end

  timeshape=size(nc.(name).data);
  if any(nc.(name).data(:)>1.E36)
   error(['There is a problem with the time data in ',name])
  end
  nc.(name).data=ridgepack_timeconvert(nc.(name).data(:),nc.(name).units,nc.(name).calendar,0);
  nc.(name).data=reshape(nc.(name).data,timeshape);

 % read WRF time data in character strings
 elseif strcmpi(name,'time') & ...
        isfield(nc,'attributes'),...
        isfield(nc.attributes,'TITLE') & ...
        ~isempty(strfind(nc.attributes.TITLE,'OUTPUT FROM WRF'))

  if(debug); disp('Extracting WRF time as a character variable'); end

  % extract the Times character string
  try
   tdimid=netcdf.inqDimID(ncid,name);
   [tdimname,tdlength] = netcdf.inqDim(ncid,tdimid);
   dimid=netcdf.inqDimID(ncid,'DateStrLen');
   [dimname,dlength] = netcdf.inqDim(ncid,dimid);
   varid=netcdf.inqVarID(ncid,'Times');
  catch
    error('Accessing time from WRF file')
  end

  % turn WRF time strings into serial time
  try
   chartime=netcdf.getVar(ncid,varid)';
   if ~isfield(nc.(name),'data')
    nc.(name).data=[1:tdlength];
   elseif (debug); 
    disp(['Converting time character record(s) ',num2str(min(nc.(name).data)),...
          ' to ',num2str(max(nc.(name).data))]);
   end
   for it=nc.(name).data
    idxc=strfind(chartime(it,1:dlength),'_');
    if ~isempty(idxc); chartime(it,idxc)=' '; end
    nc.(name).data(find(nc.(name).data==it))=datenum(chartime(it,:));
   end
   nc.(name).type='NC_DOUBLE';
  catch
   error('Converting WRF time to serial time')
  end

 % can't read it at all
 elseif strcmpi(name,'time') 

  disp(['Time units not recognized']);
  nc.(name).calendar='none'; 

 end

 % change dimension name from oldtimename to 'time' if needed
 if (isfield(nc.(name),'dimension') & any(strcmp(nc.(name).dimension,oldtimename))) & ...
     ~strcmp('time',oldtimename)
  nc.(name).dimension{strcmp(nc.(name).dimension,oldtimename)}='time';
 end

 % replace the record dimension with the label 'time' 
 if strcmp(name,oldtimename) & ~strcmp(name,'time') ; 
   nc.time=nc.(name); nc=rmfield(nc,name); 
 end

end



% ....................................................................
% Check for global attributes and fill title if attributes don't exist
if ~isfield(nc,'attributes')
 nc.attributes.title=['Data from ',ncfile];
end


% ....................................................................
% condition nc structure further if dealing with a specific model files

% Hibler-Roberts Ice-Tide model output
if isfield(nc.attributes,'title') && ...
   (~isempty(strfind(nc.attributes.title,'IARC sea ice mechanics lab : u')) | ...
    ~isempty(strfind(nc.attributes.title,'IARC sea ice mechanics lab : v'))) 
	nc=ridgepack_icetideoutput(nc);

% WRF output
elseif isfield(nc.attributes,'TITLE') & ~isempty(strfind(nc.attributes.TITLE,'OUTPUT FROM WRF'))
	nc=ridgepack_wrfoutput(nc);

% POP History file
elseif isfield(nc.attributes,'conventions') & strcmp(nc.attributes.conventions,'POP HIST conventions')
	nc=ridgepack_popoutput(nc);

% CICE History file
elseif isfield(nc.attributes,'source') & ~isempty(strfind(char(nc.attributes.source),'CICE'))
	nc=ridgepack_ciceoutput(nc);

% CPL couput output
elseif isfield(nc.attributes,'file_version') & ~isempty(strfind(char(nc.attributes.file_version),'cpl7'))
	nc=ridgepack_cploutput(nc);

end


% ....................................................................
% remove time stamp character variables written in netCDF files by ridgepack_write (if it exists)
tosstchar=0;
[variablenames,numbervariables]=ridgepack_name(nc);
for i=1:numbervariables
 name=char(variablenames{i});
 if strcmpi(name(max(1,end-6):end),'_string') & strcmp(nc.(name).type,'NC_CHAR') & ...
    isfield(nc.(name),'dimension') & any(strcmp(nc.(name).dimension,'tchar_len')) & ...
    isfield(nc.(name),'units') & strcmp(nc.(name).units,'dd-mmm-yyyy HH:MM:SS') & ...
    isfield(nc.(name),'long_name') & ~isempty(strfind(nc.(name).long_name,'values as a string'))
  nc=rmfield(nc,name);
  tosstchar=1;
 elseif (isfield(nc.(name),'dimension') & any(strcmp('tchar_len',nc.(name).dimension))) & ...
    ~strcmp('tchar_len',name)
  tosstchar=0;
 end
end
if tosstchar; 
 nc=rmfield(nc,'tchar_len'); 
 disp('Removed time units as character strings')
end

% ....................................................................
% add filename as title to the structure if none exists
if ~isfield(nc.attributes,'title')
 nc.attributes.title=ncfile(1:length(ncfile)-3);
end

% ....................................................................
% Preprocess geographical coordinates 
[dimc,coor,dims,vars]=ridgepack_content(nc);
name=intersect(vars,newvar);
if length(name)==1
 nc=ridgepack_coords(nc,char(name{1}));
else
 nc=ridgepack_coords(nc);
end
if isfield(nc,'longitude') && isfield(nc,'latitude')
 nc=ridgepack_sph2gen(nc);
end


% ....................................................................
% final conditioning of the structure
[nc,out]=ridgepack_struct(nc);

% Close the netCDF file
disp(['Closing ',ncfile]);
netcdf.close(ncid);

if debug; disp(['...Leaving ',mfilename]); end

