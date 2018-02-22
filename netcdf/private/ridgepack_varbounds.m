function [dimremove,template,tempsize,endsize,ant,time]=ridgepack_varbounds(nc,dimnames,bounds,name)

% ridgepack_varbounds - Finds reduction bounds for a specific variable in a netcdf structure
%
% function [dimremove,template,tempsize,endsize,ant,time]=ridgepack_varbounds(nc,dimnames,bounds,name)
%
% This function takes the requested bounds for reducing dimensions listed in dimnames
% and applies this list to individual variables in an nc structure that may not have all
% the dimensions listed in dimnames.  The bounds can either be in text form or in index
% form. 
%
% The main purpose of this function is to be called by the function ridgepack_reduce for
% reducing datasets. See ridgepack_reduce for further explanation.
%
% Input:
% nc       - netcdf structure (see ridgepack_struct for more details)
% dimnames - cell array of dimension variables in nc structure to be removed with 
%            slices or averages across these dimensions.
% bounds   - cell array of either upper and lower index bounds for each dimension
%            or upper and lower numeric values for each dimension over which a slice
%            or mean is to be taken.
% name     - name of variable for which these boundaries should be applied.  The
%            named variable may not possess all dimensions in dimnames, and this 
%            is taken into account by this function.
%
% Output:
% dimremove - cell array of dimensions to be removed from the variable name.
% template  - A mask template of the same dimensional size as the dimensions to be
%             removed from name.  The template is populated with NaNs for values
%             of the dimensions not to be included in slices or means, and 1s
%             for values that are to be included.
% tempsize  - This is simply size(template).
% endsize   - This is the size of of name once the dimensions in dimremove
%             have been removed from the variable as either slices or means.
% ant       - Cell array of text strings providing the boundaries over which 
%             means or slices were taken
% time      - Give the time of the slice, if a time slice is taken, otherwise
%             this remains empty.
%
% Example:
% We have the nc structure given by:
%
% netcdf structure nc {
%  attributes {global}
%     Conventions: 'CF-1.0'
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
%  time [01-Jan-1996 00:00:00 to 30-Jun-1997 18:00:00, mean step 06:00:00.000]
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
% We want to reduce the dimensions by removing time, taking a time slice at 01-Jan-1996 12:00:00.
% The dimension to be removed is 'time', so dimnames={'time'}. Since we want to take a time slice
% bounds={{'01-Jan-1996 12:00:00','01-Jan-1996 12:00:00'}}, although, knowing that the timestep
% is 6 hours, we could also have entered  bounds={[3 3]} since this slice is the third value
% in the time vector.   If we had wanted to take a mean over, say, the second day of the record
% instead, we would have entered either bounds={{'02-Jan-1996 00:00:00','02-Jan-1996 18:00:00'}}
% or bounds={[5 8]}.  The advantage of the former value is that we do not need to know the 
% timestep of the data, so we could have entered  bounds={{'01-Jan-1996 22:00:00','02-Jan-1996 20:00:00'}}
% to obtain the same result, since this still provides the upper and lower bounds for the 
% data at 6-hourly resolution.
%
% Now, we want to know what the bounds are for msl, so name='msl'. The output for the case of 
% bounds={[3 3]} would be the following:
%
% dimremove={'time'}
% template=[nan; nan; 1; nan; ...; nan] (this vector has 2188 elements, all NaNs except the third)
% tempsize=size(template)=[2188 1]
% endsize=[73 144] (the size of the reduced msl variable when the bounds are applied to it).
% ant={'01-Jan-1996 12:00:00'}
%
% Note that if we were to set name='latitude' or name='longitude' with all other settings
% remaining the same as above, we would get the following output:
%
% dimremove={}
% template=nan
% tempsize=[1 1]
% endsize=[] 
% ant={}
%
% The reason being that 'time' is not a dimension of either latitude or longitude, and
% so there would be no bounds applied to them. 
%
% Andrew Roberts, Naval Postgraduate School, March 2018  (afrobert@nps.edu)
%

global debug;
if debug; disp(['Entering ',mfilename,'...']); end

% default
time=[];

% Dimensions to be removed from this variable
m=0;
dimremove={};
for j=1:length(dimnames)
	if any(strcmp(char(dimnames{j}),nc.(name).dimension))
		m=m+1;
		dimremove{m}=dimnames{j};
	end
end

% get the length of the dimensions descriptor
dimn=length(dimremove);

% get the length of the actual data dimensions to be removed that
% takes into account trailing singleton dimensions being removed by permute
dimred=max(0,length(dimremove)-(length(nc.(name).dimension)-ndims(nc.(name).data)));

% Check that dimensions are being removed for var, else simply
% write to the new structure.
if isempty(dimremove)
         if debug; disp(['No dimensions being removed from ',name]); end
	 ncr.(name)=nc.(name);
	 template=NaN;
	 tempsize=size(squeeze(template));
 	 endsize=[];
	 ant={};
	 return;
end

% Obtain shape of the array being reduced
startsize=size(nc.(name).data);
if debug; disp([name,' initial shape is ',num2str(startsize)]); end

% Check that not all dimensions are being removed
if isempty(setdiff(nc.(name).dimension,dimremove))
         disp(['Removing all dimensions of ',name]);
elseif length(startsize)<=dimn & length(nc.(char(nc.(name).dimension(end))).data)>1
	 error(['The number of dimensions being removed for ',name,...
	        ' must be less than its total dimensions']);
elseif dimn>4
         error('The number of dimnsions being removed must not exceed 4');
end


% Calculate shape of final data array
endsize=[startsize(1:end-dimred) 1];
if debug; disp([name,' final shape is ',num2str(endsize)]); end

% Calculate the shape of the template and fill with NaNs
if dimred==0
	tempsize=[1 1];
else
	tempsize=[startsize(end-dimred+1:end) 1];
end

% This line wrong - replaced with next one below
%template=squeeze(ones(tempsize)); 
template=ones(tempsize);

% Prepare dummy bounds array containing only the boundaries 
% for the dimensions being removed from this particular variable
% Also check the array bounds when they have beeb explictly provided
% and prepare the long_name text for each dimension.

bo={};
cond={};
m=0;
for j=1:length(dimnames)
         
         % check bounds are not a mixture of text and indices
	 if ((isnumeric(bounds{j}(1)) & ~isnumeric(bounds{j}(2))) | ...
	    (isnumeric(bounds{j}(2)) & ~isnumeric(bounds{j}(1))))
	    error('boundaries are a mixture of text and indices')
         end

	 % if this dimension is to be removed, add value to numeric bo index array
         if any(strcmp(char(dimnames{j}),dimremove))

	  m=m+1; % increment the index of dimensions being removed from this variable

	  % if bounds are numeric, simply add value to numeric array, checking value is in range
	  if isnumeric(bounds{j})

	   bo{m}=bounds{j};
	   
	   % check values are in range and reset, or stop if dimension is time to prevent runaway
           if bo{m}(1)<1
	    disp(['First index ',num2str(bo{m}(1)),' for ',char(dimnames{j}),' too small']);
	    bo{m}(1)=1;
	    if strcmp('time',char(dimnames{j}))
   	     error('Check time indices');
            else
             disp(['Reset to minumum possible value ',num2str(bo{m}(1))]);
            end
           end
           if bo{m}(1)>length(nc.(char(dimnames{j})).data)
	    disp(['First index ',num2str(bo{m}(1)),' for ',char(dimnames{j}),' too large']);
	    bo{m}(1)=length(nc.(char(dimnames{j})).data);
	    if strcmp('time',char(dimnames{j}))
   	     error('Check time indices');
            else
             disp(['Reset to maximum possible value ',num2str(bo{m}(1))]);
            end
           end
           if bo{m}(2)<1
	    disp(['Second index ',num2str(bo{m}(2)),' for ',char(dimnames{j}),' too small']);
	    bo{m}(2)=1;
	    if strcmp('time',char(dimnames{j}))
   	     error('Check time indices');
            else
             disp(['Reset to minumum possible value ',num2str(bo{m}(2))]);
            end
           end
           if bo{m}(2)>length(nc.(char(dimnames{j})).data)
	    disp(['Second index ',num2str(bo{m}(2)),' for ',char(dimnames{j}),' too large']);
	    bo{m}(2)=length(nc.(char(dimnames{j})).data);
	    if strcmp('time',char(dimnames{j}))
   	     error('Check time indices');
            else
             disp(['Reset to maximum possible value ',num2str(bo{m}(2))]);
            end
           end
           if bo{m}(1)>bo{m}(2)
            a=bo{m}(2);
            bo{m}(2)=bo{m}(1);
	    bo{m}(1)=a;
	    disp(['Switched index values for ',char(dimnames{j})]);
           end

	  % translate numeric boundaries on dimension to indices
          else

           if strcmp('time',char(dimnames{j}))
            lower=datenum(char(bounds{j}(1)));
            upper=datenum(char(bounds{j}(2)));
           else
            lower=str2num(char(bounds{j}(1)));
            upper=str2num(char(bounds{j}(2)));
           end
	   if isempty(lower) | isempty(upper)
            error(['Boundaries input for ',char(dimnames{j}),' are incorrect']);
           end
	   if lower>upper
	    a=upper;
	    upper=lower;
	    lower=a;
	   end
	   cond{m}=['nc.',char(dimnames{j}),'.data>=',num2str(lower),...
	            ' & nc.',char(dimnames{j}),'.data<=',num2str(upper)];
           disp(['Looking for values ',char(cond{m})]);
           eval(['a=',char(cond{m}),';']);

	   % check requested numerical boundaries contains data
	   if ~any(a)
            disp(['Specification ',char(cond{m}),' contains no data']);
	    val=mean(upper,lower);
	    bo{m}(1) = find(abs(nc.(char(dimnames{j})).data(:)-val)==...
	                min(abs(nc.(char(dimnames{j})).data(:)-val)));
            bo{m}(2)=bo{m}(1);
	    if strcmp('time',char(dimnames{j}))
   	     error('Check time indices');
            else
             disp(['Reset to index ',num2str(bo{m}(1)),' for which ',char(dimnames{j}),...
	           ' = ',num2str(nc.(char(dimnames{j})).data(bo{m}(1)))]);
            end
           else
            bo{m}=[min(find(a)) max(find(a))];
	   end 

	  end 
	 end 
end 

% check that bo has been correctly filled
if isempty(bo);
        error('Unable to allocate values to boundaries array bo');
end

% prepare text for output in metadata of slice/means taken
ant=cell(length(dimremove),1);
time=[];
for j=1:length(dimremove)

	 if strcmp('time',char(dimnames{j})); 
	  if bo{j}(1)==bo{j}(2)
	   ant{j}=[datestr(nc.(char(dimremove{j})).data(bo{j}(1)),0)];
	   time=nc.(char(dimremove{j})).data(bo{j}(1));
	  else
           if nc.(char(dimremove{j})).data(bo{j}(2))>...
              nc.(char(dimremove{j})).data(bo{j}(1))
	    ant{j}=['[',datestr(nc.(char(dimremove{j})).data(bo{j}(1)),0),'>',...
	            datestr(nc.(char(dimremove{j})).data(bo{j}(2)),0),']'];
	   else
	    ant{j}=['[',datestr(nc.(char(dimremove{j})).data(bo{j}(2)),0),'>',...
	            datestr(nc.(char(dimremove{j})).data(bo{j}(1)),0),']'];
	   end
	   time=mean(nc.(char(dimremove{j})).data(bo{j}(1:2)));
	  end
	 else
	  units=ridgepack_units(nc,char(dimnames{j})); % output units for plotting
	  if ~isfield(nc.(char(dimnames{j})),'units') | ...
	      isempty(strfind(nc.(char(dimnames{j})).units,'degrees')) 
           param=[char(dimnames{j}),'='];
          else
           param='';
   	  end
	  if bo{j}(1)==bo{j}(2)
	   ant{j}=[param,num2str(nc.(char(dimnames{j})).data(bo{j}(1))),units];
	  else
           if nc.(char(dimremove{j})).data(bo{j}(2))>...
              nc.(char(dimremove{j})).data(bo{j}(1))
	    ant{j}=['[',param,num2str(nc.(char(dimremove{j})).data(bo{j}(1))),':',...
	            num2str(nc.(char(dimremove{j})).data(bo{j}(2))),units,']'];
	   else
	    ant{j}=['[',param,num2str(nc.(char(dimremove{j})).data(bo{j}(2))),':',...
	            num2str(nc.(char(dimremove{j})).data(bo{j}(1))),units,']'];
	   end
	  end
	 end
end


% write boundaries to screen for information purposes only
% run checks on boundaries at the same time
char0=(['Reducing dimensions for ',name]);
if dimn==1
	  char1=[' along indices ',char(dimremove{1}),'=',num2str(bo{1}(1)),':',num2str(bo{1}(2))];
	  disp([char0,char1]);
	  if length(bo)~=1;error('bo length should be 1');end;
elseif dimn==2
	  char1=[' along indices ',char(dimremove{1}),'=',num2str(bo{1}(1)),':',num2str(bo{1}(2))];
	  char2=[' ',char(dimremove{2}),'=',num2str(bo{2}(1)),':',num2str(bo{2}(2))];
	  disp([char0,char1,char2]);
	  if length(bo)~=2;error('bo length should be 2');end;
elseif dimn==3
	  char1=[' along indices ',char(dimremove{1}),'=',num2str(bo{1}(1)),':',num2str(bo{1}(2))];
	  char2=[' ',char(dimremove{2}),'=',num2str(bo{2}(1)),':',num2str(bo{2}(2))];
	  char3=[' ',char(dimremove{3}),'=',num2str(bo{3}(1)),':',num2str(bo{3}(2))];
	  disp([char0,char1,char2,char3]);
	  if length(bo)~=3;error('bo length should be 3');end;
elseif dimn==4
	  char1=[' along indices ',char(dimremove{1}),'=',num2str(bo{1}(1)),':',num2str(bo{1}(2))];
	  char2=[' ',char(dimremove{2}),'=',num2str(bo{2}(1)),':',num2str(bo{2}(2))];
	  char3=[' ',char(dimremove{3}),'=',num2str(bo{3}(1)),':',num2str(bo{3}(2))];
	  char4=[' ',char(dimremove{4}),'=',num2str(bo{4}(1)),':',num2str(bo{4}(2))];
	  disp([char0,char1,char2,char3,char4]);
	  if length(bo)~=4;error('bo length should be 4');end;
end

% Fill template with ones wherever calculations are to be done, else with NaNs 
template(:)=NaN;
if dimn==1
	  template(bo{1}(1):bo{1}(2))=1;
elseif dimn==2
	  template(bo{1}(1):bo{1}(2),bo{2}(1):bo{2}(2))=1;
elseif dimn==3
	  template(bo{1}(1):bo{1}(2),bo{2}(1):bo{2}(2),bo{3}(1):bo{3}(2))=1;
elseif dimn==4
	  template(bo{1}(1):bo{1}(2),bo{2}(1):bo{2}(2),...
	       bo{3}(1):bo{3}(2),bo{4}(1):bo{4}(2))=1;
else
	  error(['The number of dimensions to be reduced exceeds the limit of 4']);
end

% check that template is not all NaNs 
if isempty(template(~isnan(template)))
	       disp(['The template for ',name,' is blank (filled with NaNs)']);
	       error('Check the bounds are within range of the data')
end

if debug; disp(['...Leaving ',mfilename]); end

