function [nc]=ridgepack_groom(nc, name)

% ridgepack_groom - Grooms a netcdf structure for writing or after reading netcdf files
%
% function [nc]=ridgepack_groom(nc, name)
% 
% This function grooms individual components of a netcdf structure (see 
% ridgepack_struct for explanation of this structure and required inputs).  The
% fields are ordered into a logical construction for listing on screen. A
% fill value is added if NaN's are present but no fill value is 
% present.  Missing values are transferred to fill values and removed 
% from the structure, as is the scale factor and offset fields not used 
% for writing or within matlab except on reading the original netcdf file
% by nctools. The type and units of the data are checked.  Time units 
% are kept in place for writing a netcdf file, and data type is that 
% of a netcdf write, but neither are representative of the data in matlab. 
%
% If the type of the data is nc_char or if a grid mapping is set, the 
% normal checks on the data are bypassed.
%
% INPUT:
%
% nc   - netcdf data structure input
% name - name of component being groomed in the netcdf structure.
%
%
% OUTPUT:
%
% nc   - groomed data structure for nc.name
%
% Ridgepack Version 1.0
% Andrew Roberts, Naval Postgraduate School, March 2018 (afrobert@nps.edu)
%

global debug;
if debug; disp(['Entering ',mfilename,'...']); end

% Check data is a structure
if not(isstruct(nc));
    error([inputname(1),' is not a structure']);
end     

% set conditions for which a bypass of normal data checks is possible
if isfield(nc.(name),'type') && strcmp(nc.(name).type,'NC_CHAR');
        bypass=true;
elseif isfield(nc.(name),'grid_mapping_name') 
        disp(['Changing ',name,' type to NC_CHAR']);
        bypass=true;
else
	bypass=false;
end

if ~isfield(nc.(name),'type') & isfield(nc.(name),'data')
  nc.(name).type=ridgepack_invert(class(nc.(name).data));
end

% check that a long name has been provided, and clean it up for using in
% plotting by removing inadvertent TeX descriptors
if isfield(nc.(name),'long_name') ; 
	if not(isstrprop(nc.(name).long_name, 'alphanum'));
	 disp(['long_name is not alphanumeric for ',name,':',nc.(name).long_name]);
	 nc.(name).long_name=name;
	else
	 pos_=findstr(nc.(name).long_name,'_');
	 if ~isempty(pos_)
	  for i=1:length(pos_)
	   nc.(name).long_name(pos_(i))=' ';
	  end
	 end
	end
elseif isfield(nc.(name),'description') ;
	nc.(name).long_name=nc.(name).description;
	nc.(name)=rmfield(nc.(name),'description');
elseif ~bypass
	if ~isstruct(nc.(name)); 
         disp(['nc.',name,' appears to be an incorrectly assigned variable.']);
	 disp(['It appears that the nc structure has been corrupted']);
	 error(['nc.',name,' is not a structure.'])
        end
	nc.(name).long_name=name ; 
end

% reset time units
if (strcmpi(name,'time') | strcmpi(name,'time_bounds')) & ...
   isfield(nc.(name),'type') & ~strcmp(nc.(name).type,'NC_CHAR') 

        % find NaN records in time 
        if isfield(nc.(name),'data')
         removeindex=find(isnan(nc.(name).data(:)));
        else
         error('Unable to find any time data')
        end

	% get minumum time
        try
	 tvec=datevec(min(nc.(name).data(:)));
        catch
	 if debug; disp(datestr(nc.(name).data)); end
         warning(['Time specification error. Minimum time=',num2str(min(nc.(name).data(:)))])
         error('Time specification')
        end 

	if debug; disp(['Date minimum is ',datestr(datenum(tvec))]); end

        % set a calendar to gregorian if none is set
        if ~isfield(nc.(name),'calendar'); 
	 nc.(name).calendar='proleptic_gregorian'; 
	end

        % check calendar status
        if any(strcmpi(nc.(name).calendar,{'standard','gregorian','noleap','365_day'}))
         caltype=1;
        elseif strcmpi(nc.(name).calendar,'proleptic_gregorian')
         caltype=2;
        elseif any(strcmpi(nc.(name).calendar,{'all_leap','366_day','julian','360_day'}))
         caltype=3;
	else
         caltype=0;
        end

	% set default time coordinates here for conversion when writing
	% to netcdf from the serial time represented in matlab
        if isempty(tvec)
         disp('You are missing time data')
         error('The time vector is empty')
	elseif tvec(1)>=2000 & (caltype==1 | caltype==2)
 	 nc.(name).units='hours since 2000-01-01 00:00:0.0';
	elseif tvec(1)>=1900 & (caltype==1 | caltype==2)
 	 nc.(name).units='hours since 1900-01-01 00:00:0.0';
	elseif tvec(1)>=1800 & (caltype==1 | caltype==2)
 	 nc.(name).units='hours since 1800-01-01 00:00:0.0';
	elseif tvec(1)>=1700 & (caltype==1 | caltype==2)
 	 nc.(name).units='hours since 1700-01-01 00:00:0.0';
	elseif tvec(1)>=1600 & (caltype==1 | caltype==2)
 	 nc.(name).units='hours since 1600-01-01 00:00:0.0';
	elseif datenum(tvec)>=datenum([1582 10 15]) & caltype==1
 	 nc.(name).units='hours since 1582-10-15 00:00:0.0';
	elseif caltype==2 | caltype==3
 	 nc.(name).units=['days since ',datestr(fix(datenum(tvec)),'yyyy-mm-dd')];
	end

end

% make sure fill values have been assigned to latitude and longitude
if strcmp(name,'latitude')
        nc.(name).data(nc.(name).data>90)=NaN;
        nc.(name).data(nc.(name).data<-90)=NaN;
end

% make sure fill values have been assigned to latitude and longitude
if strcmp(name,'longitude')
        nc.(name).data(nc.longitude.data>360)=NaN;
        nc.(name).data(nc.longitude.data<-360)=NaN;
        nc.(name).data(isnan(nc.latitude.data))=NaN;
end


% remove units if it is empty
if isfield(nc.(name),'units') & isempty(nc.(name).units) ;
	nc.(name)=rmfield(nc.(name),'units');
end


% check that dimensions are provided where needed
if isfield(nc.(name),'dimension') & ~isempty(nc.(name).dimension); 
	if not(iscell(nc.(name).dimension)); 
 	 error(['Dimension info not a cell array for ',name]);
	end
elseif ~strcmp(nc.(name).type,'NC_CHAR') & isfield(nc.(name),'data') & ...
        length(nc.(name).data)>1
	disp(['No dimension provided for ',name]) ; 
	nc.(name).dimension={};
end

% check the data
if isfield(nc.(name),'data') & ~isempty(nc.(name).data); 
	if not(isempty(nc.(name).data(isinf(nc.(name).data))));
	 disp(['WARNING: There are infinites in the data of ',name]);
	 disp(['         Converting infinites to NaNs']);
         nc.(name).data(isinf(nc.(name).data))=NaN;
	end
        if isempty(nc.(name).data) ;
	 disp(['WARNING: The data array is empty for ',name]) ;
        end
elseif ~bypass
	disp(['WARNING: No data array provided for ',name]) ;
        nc.(name).data=[];
end

% move nc.(name).FillValue to nc.(name).fillvalue
if isfield(nc.(name),'FillValue'); 
	nc.(name).fillvalue=nc.(name).FillValue;
	nc.(name)=rmfield(nc.(name),'FillValue');
end

% check for missing values in the data for discrepencies
if isfield(nc.(name),'fillvalue'); 

	if not(isnumeric(nc.(name).fillvalue)) ;
	 error(['The fill value in ',name,' is non-numeric']);
        end

	if ~any(isnan(nc.(name).data(:)));
         nc.(name)=rmfield(nc.(name),'fillvalue');
        end

elseif ~bypass & any(isnan(nc.(name).data(:)));

	nc.(name).fillvalue=min(-99999,floor(9*min(nc.(name).data(:))));

	% This line removed to stop filling fill values with NaN unless being written
	%nc.(name).data(isnan(nc.(name).data))=nc.(name).fillvalue;

	if debug
	 disp(['WARNING: There were NaNs but no assigned fill value for ',name]);
	 disp(['         The fill value ',num2str(nc.(name).fillvalue),' was assigned.']);
	end

	% make sure the class of the fill value is identical to the data
	mclass=ridgepack_onvert(nc.(name).type);
	eval(['nc.',name,'.fillvalue=',mclass,'(nc.',name,'.fillvalue);'])

end


% remove fields that are non-standard in the CF convention, but frequently occur
if isfield(nc.(name),'missing_value') ;
	nc.(name)=rmfield(nc.(name),'missing_value');
end

if isfield(nc.(name),'actual_range') ;
	nc.(name)=rmfield(nc.(name),'actual_range');
end

% remove range fields (these will be added if necessary when writing a netcdf file)
if isfield(nc.(name),'valid_range') ;
	nc.(name)=rmfield(nc.(name),'valid_range');
end

if isfield(nc.(name),'valid_max') ;
	nc.(name)=rmfield(nc.(name),'valid_max');
end

if isfield(nc.(name),'valid_min') ;
	nc.(name)=rmfield(nc.(name),'valid_min');
end

% remove invalid axis definitions for multiple dimension variables
if isfield(nc.(name),'axis') & isfield(nc.(name),'dimension') & length(nc.(name).dimension)>1 ;
	nc.(name)=rmfield(nc.(name),'axis');
end
 
if debug; disp(['...Leaving ',mfilename]); end

