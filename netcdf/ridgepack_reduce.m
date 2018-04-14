function [ncr,nc]=ridgepack_reduce(nc,dimnames,bounds,varmean)

% ridgepack_reduce - Reduce dimensions in an nc structure by slicing or taking weighted means 
%
% function [ncr,nc]=ridgepack_reduce(nc,dimnames,bounds,varmean)
%
% This is a sophisticated function that removes dimensions from an nc structure and 
% recognizes calculations that need to be performed to do this for geolocated and masked
% gridded data without user prompting. Operations performed during execution are verbose
% so the user understands exactly what automated operations are taking place.
%
% INPUT:
%
%
% nc        - nc structure (see the help page of ridgepack_struct for more information)
%
% dimnames  - cell array of dimension variable names and order of dimensions to be 
%             removed from var. These may be specified in any order, with ridgepack_reduce 
%             being clever enough to work out the way in which they will be removed 
%             from each variable. So for a variable with dimensions 
%             {'time'  'latitude'  'longitude'} for which you want a latitudinal
%             average, you would set dimnames={'time','longitude'} but 
%             dimnames={'longitude','time'} would work in exactly the same way.
%
% bounds    - cell array of bounds to be applied to the dimension in reducing the data.  
%             This input may take several forms:
%             1) The boundaries may be specified in terms of index values:
%             e.g. dimnames={'time'}, bounds={[1 4]} will take a linearly 
%             weighted-in-time mean between time slices 1 and 4 in the nc structure.
%             2) The boundaries may be specified in terms of dimension values:
%             e.g. dimnames={'time'}, bounds={{'01-Jan-2007','04-Jan-2007 06:00:00.0'}}
%             which will take a linearly weighted-in-time mean between 
%             Jan 1 12AM  and Jan 4 6AM 2007
%             3) Multiple combinations where multiple dimensions are to be removed 
%             at once: e.g. dimnames={'time','latitude'}, bounds={[1 4],{'70 90'}} 
%             for latitude in degrees will take a linearly weighted-in-time mean 
%             and an area weighted-on-a-sphere meanfor grid points between 70 and 90 
%             North. Note that, as shown in the last example, that the bounds for 
%             reducing a dimension appears in the same order as the name of the 
%             dimensions in dimnames. There must be two elements in bounds, even if 
%             they are the same to simply take a slice.
%
% varmean   - logical set to 'true' to gain more statistical output rather than 
%             just means during the reduction process. The default is 'false'.  
%             This is an optional argument.
%
% CICEpatch - Patch CICE monthly mean data that is missing a time_bounds descriptor
%
%
%
% OUTPUT:
%
% ncr  - Reduced netcdf structure with requested dimensions removed and variable 
%        means, and, if varmean is true, the standard deviations and number of 
%        samples in means, and effective sample size.  These are calculated based
%        on the weightings of each sample, which, if done in time, depends in the 
%        time_bounds or difference in time between samples if the times are less than
%        28 days apart.
%
% nc   - The original nc structure with added information during computations that can be 
%        reused in another reduction operation to reduce future computation times.
%
% Implicit functionality of ridgepack_reduce:
%
% * In addition simply taking means or slides through dimension of data in nc structures, 
% this function recognizes where area-weighted means need to be taken, both for straight 
% spherical grids, and for more complicated limited area or generalised curvilinear grids
% by taking a line integral around each grid point on the sphere to calculate the area
% within each grid point.   Note that it is assumed that grid cells have four vertices 
% or less. Since the line integral method takes time, it is advisable to retain the 
% information from this calculation done the first time which is added to nc in the 
% output. Taking the line integrals on a geoid requires a simple one-line change in 
% the code, but this is not provided in this code.
%
% * ridgepack_reduce searches the nc structure for a variable called "mask" and uses this 
% to mask any means being taken, where it is assumad the mask is binary, with 0 being 
% the points on a grid not to be included in calculations, and 1 being points that are.
%
% * You may specify only the dimension name if the dimension being removed only has 
% a single value and you are only reducing the structure by one dimension
%
% Example:
%
% Here is an example of reading in a netcdf file, and reducing the structure of an 
% ECMWF ERA40 file to remove time and longitude dimensions using 
% dimnames={'longitude','time'} and bounds={[1 144],{'01-Jan-1958','03-Jan-1958'}}. 
% The structure is first imported from a netcdf file called era40.1958.nc, using 
% nc=ridgepack_clone('era40.1958.nc','p2t'), and then is reduced using 
% [ncr,nc]=ridgepack_reduce(nc,{'longitude','time'},{[1 144],{'01-Jan-1958','03-Jan-1958'}}):
%
% >> nc=ridgepack_clone('era40.1958.nc','p2t')
%
% Gives the initial nc structure:
%  
% netcdf structure nc {
%  attributes {global}
%     Conventions: 'CF-1.0'
%         history: '2007-11-12 20:00:39 GMT by mars2netcdf-0.92'
%           title: 'era40.1958'
%  
%  4 variables, dimension=[ ], data=( )
%  
%  latitude [-90 to 90, range 180, step -2.5]
%          data: [73x1 double]
%     dimension: {'latitude'}
%     long_name: 'latitude'
%          type: 'NC_FLOAT'
%         units: 'degrees_north'
%  
%  longitude [0 to 357.5, range 357.5, step 2.5]
%          data: [144x1 double]
%     dimension: {'longitude'}
%     long_name: 'longitude'
%          type: 'NC_FLOAT'
%         units: 'degrees_east'
%  
%  time [01-Jan-1958 00:00:00 to 01-Jan-1959 00:00:00, timestep 06:00:00.000]
%          data: [1461x1 double]
%     dimension: {'time'}
%     long_name: 'time'
%          type: 'NC_FLOAT'
%         units: 'hours since 1900-01-01 00:00:0.0'
%  
%  p2t (193.1944 to 321.9842, range 128.7898, median 282.9082)
%     coordinates: 'latitude longitude'
%            data: [1461x73x144 double]
%       dimension: {'time'  'latitude'  'longitude'}
%       fillvalue: -99999
%       long_name: '2 metre temperature'
%            type: 'NC_FLOAT'
%           units: 'K'
%  
% }
% 
%  
% >> [ncr,nc]=ridgepack_reduce(nc,{'longitude','time'},{[1 144],{'01-Jan-1958','03-Jan-1958'}});
%
% Gives the reduced nc structure ncr and the initial nc structure with the area of 
% cells included:
% 
% netcdf structure ncr {
%  attributes {global}
%     Conventions: 'CF-1.0'
%         history: '2007-11-12 20:00:39 GMT by mars2netcdf-0.92'
%           title: 'era40.1958 reduced dataset'
%         comment: 'Reduced dataset via slices or means through: longitude time'
%  
%  3 variables, dimension=[ ], data=( )
%  
%  latitude [-90 to 90, range 180, step -2.5]
%          data: [73x1 double]
%     dimension: {'latitude'}
%     long_name: 'latitude'
%          type: 'NC_FLOAT'
%         units: 'degrees_north'
%  
%  p2t (234.9903 to 299.8144, range 64.824, delta -2.5)
%          data: [73x1 double]
%     dimension: {'latitude'}
%     long_name: '2 metre temperature MEAN [0:357.5^{o}E] [01-Jan-1958 00:00:00>02-Jan-1958 18:00:00] N=1152'
%          type: 'NC_FLOAT'
%         units: 'K'
%  
%  p2t_std (0.13339 to 12.8638, range 12.7304, delta -2.5)
%          data: [73x1 double]
%     dimension: {'latitude'}
%     long_name: '2 metre temperature STDEV [0:357.5^{o}E] [01-Jan-1958 00:00:00>02-Jan-1958 18:00:00] N=1152'
%          type: 'NC_FLOAT'
%         units: 'K'
%
%  ...
%  ...
%  ...
%  
% }
%  
% netcdf structure nc {
%  attributes {global}
%     Conventions: 'CF-1.0'
%         history: '2007-11-12 20:00:39 GMT by mars2netcdf-0.92'
%           title: 'era40.1958'
%  
%  5 variables, dimension=[ ], data=( )
%  
%  latitude [-90 to 90, range 180, step -2.5]
%          data: [73x1 double]
%     dimension: {'latitude'}
%     long_name: 'latitude'
%          type: 'NC_FLOAT'
%         units: 'degrees_north'
%  
%  longitude [0 to 357.5, range 357.5, step 2.5]
%          data: [1x144 double]
%     dimension: {'longitude'}
%     long_name: 'longitude'
%          type: 'NC_FLOAT'
%         units: 'degrees_east'
%  
%  time [01-Jan-1958 00:00:00 to 01-Jan-1959 00:00:00, timestep 06:00:00.000]
%          data: [1461x1 double]
%     dimension: {'time'}
%     long_name: 'time'
%          type: 'NC_FLOAT'
%         units: 'hours since 1900-01-01 00:00:0.0'
%        weight: [1x1461 double]
%  
%  cell_area (8.2564e-07 to 0.00015166, range 0.00015083, median 0.0001071)
%     coordinates: 'latitude longitude'
%            data: [73x144 double]
%       dimension: {'latitude'  'longitude'}
%       long_name: 'fractional area on the sphere'
%            type: 'NC_DOUBLE'
%  
%  p2t (193.1944 to 321.9842, range 128.7898, median 282.9082)
%     coordinates: 'latitude longitude'
%            data: [73x144x1461 double]
%       dimension: {'latitude'  'longitude'  'time'}
%       fillvalue: -99999
%       long_name: '2 metre temperature'
%            type: 'NC_FLOAT'
%           units: 'K'
%  
% }
%            
%
% Ridgepack Version 1.0
% Andrew Roberts, Naval Postgraduate School, March 2018 (afrobert@nps.edu)

global debug;
%debug=true;

if debug; disp(['Entering ',mfilename,'...']); end

% Check data is a structure
if not(isstruct(nc));
 error([inputname(1),' is not a structure']);
end

% Check that dimensions are specified
if nargin<2
 error('no dimensions specified')
elseif ~iscell(dimnames)
 error('dimensions must be in a cell array');
elseif isempty(dimnames)
 if debug; disp('Dimensions cell array is empty, nc structure unchanged'); end
 ncr=nc;
 return
end

% If removing singleton dimensions, no need to specify bounds
if nargin<3 | (iscell(bounds) & isempty(bounds))
 if debug; disp('Reducing over the entire dimensional arrays requested'); end
 bounds=cell(1,length(dimnames));
 for i=1:length(dimnames)
  if min(size(nc.(char(dimnames{i})).data))>1 | ...
     length(size(nc.(char(dimnames{i})).data))>2
   error(['The bounds of ',char(dimnames{i}),' suggest it is not a dimension variable.'])
  end
  bounds{i}=[1 max(size(nc.(char(dimnames{i})).data))];
  if debug;
   disp([char(dimnames{i}),' bounds are [1 ',...
     num2str(max(size(nc.(char(dimnames{i})).data))),']']); 
  end
 end
end

% Turn off standard deviation calculations if not requested
if nargin<4
 varmean=false;
 meanonly=false;
elseif ~strcmp(class(varmean),'logical')
 error('varmean must be a logical')
end

% Check input
if  ~iscell(bounds)
 error('bounds not supplied in cell array');
elseif length(dimnames)~=length(bounds)
 error('bounds must have upper and lower bounds for each element of dimnames')
elseif isempty(dimnames)
 disp('No reduction performed to the netcdf structure');
 ncr=nc;
 return
end

% check bounds to see if we can fastrack...
fasttrack=true;
for i=1:length(bounds)
  if length(bounds{i})~=2 
    error(['Incorrect number of elements in bounds ',num2str(i)]);
  end
  if fasttrack 
   if ~isnumeric(bounds{i}(1)) & ~isnumeric(bounds{i}(1)) & ...
      ~strcmp(char(bounds{i}(1)),char(bounds{i}(2)))
    fasttrack=false;
   elseif isnumeric(bounds{i}(1)) & isnumeric(bounds{i}(1)) & ...
      ~(bounds{i}(1)==bounds{i}(2))
    fasttrack=false;
   end
  end
end

if debug
 if fasttrack
  disp('Taking a single slice of data...')
 else
  disp('Taking more than just a slice of cake...');
 end
end

% make sure coordinate have been assigned
nc=ridgepack_coords(nc);

% remove mesh variables and dimensions from the structure because
% they are of no use in reducing the dimensions of data variables
[dimc,coor,dims,vars,mdim,mesh,supp]=ridgepack_content(nc);

if isempty(ridgepack_union(mdim,ridgepack_union(dimc,dims))); ...
  error('nc contains no dimensions to be reduced'); ...
end
if isempty(ridgepack_union(mesh,ridgepack_union(coor,vars))); ...
  error('nc contains no variables to be reduced'); ...
end

% check that dimensions to be removed exist in the dataset
if isempty(intersect(dimnames,ridgepack_union(mdim,ridgepack_union(dimc,dims))));
 error(['No specified dimensions match variable dimensions: ',ridgepack_cellcat(dimnames)])
end

% Set up weighting arrays in variables
for i=1:length(vars)
	name=char(vars{i});
	nc.(name).weight=ones([size(nc.(name).data)]);
	nc.(name).weight(isnan(nc.(name).data))=0;
end

if ~fasttrack

 % Calculate the cell areas if any coordinate dimensions are being removed
 if ~isempty(intersect(coor,dimnames)) 

	nc=ridgepack_cellarea(nc);
	nc=ridgepack_shuffle(nc,nc.cell_area.dimension);

	% apply these area weigthings to each variable
	for i=1:length(vars)
 	name=char(vars{i});
	 if ~isempty(intersect(nc.cell_area.dimension,nc.(name).dimension))
		snc=size(nc.(name).data);
		numd=ndims(nc.cell_area.data);
		if size(nc.(name).weight)~=snc
	 	 error(['1: weight and data arrays different size for ',name])
		end
		newsize=[max(prod(snc(1:end-numd)),1),max(prod(snc(end-numd+1:end)),1)];
		nc.(name).weight=reshape(nc.(name).weight,newsize);
		for j=1:max(prod(snc(1:end-numd)),1)
		 nc.(name).weight(j,:)=nc.(name).weight(j,:).*nc.cell_area.data(:)';
		end
		nc.(name).weight=reshape(nc.(name).weight,snc);
		disp(['Adding area-weighted mean to ',name]);
	 end
	end

 % If cell area was not calculated but coors variables exist
 elseif ~isempty(setdiff(coor,{'latitude','longitude'})) 
	disp(['Don''t know how to calculate weightings for coordinates: ',...
             ridgepack_cellcat(setdiff(coor,{'latitude','longitude'}))]);
	%error('weigtings not calculated for coordinates provided');

 end

 % Create linear weightings for non-coordinate dimensions in dims
 % Weightings are kept in the nc structure, but not written to output.
 for i=1:length(dims)
  name=char(dims{i});

  if any(strcmpi(name,dimnames))

   if isfield(nc.(name),'weight')
	disp(['Using pre-existing weightings for ',name]);
   else

       if debug; disp(['Starting to create weightings for ',name]); end

       l=length(nc.(name).data);

       if l>1

        if strcmpi(name,'time') & isfield(nc,'time_bounds')

         if length(nc.time_bounds.dimension)>2
          error('time_bounds has more than two dimensions')
         else
          disp('Creating non-linear time weightings from time_bounds')
         end

         % convert to netcdf time units so as to take calendar type into account
         timeshape=size(nc.time_bounds.data);
         time_bounds=ridgepack_timeconvert(nc.time_bounds.data(:),...
                                   'days since 0000-01-01 00:00:0.0',...
                                   nc.time_bounds.calendar,1);
         time_bounds=reshape(time_bounds,timeshape);
         d2place=find(~strcmpi(nc.time_bounds.dimension,'time'));
         timeplace=find(strcmpi(nc.time_bounds.dimension,'time'));
         time_bounds=permute(time_bounds,[d2place timeplace]);

         % take difference to get weightings
         weight=diff(time_bounds,1,1);

        elseif strcmpi(name,'time') 

         disp('WARNING: No time bounds to calculate time weights')

         % convert to netcdf time units so as to take calendar type into account
         if length(nc.time.dimension)>1
          error('times has more than two dimensions')
         end

         % convert to netcdf time units so as to take calendar type into account
         time=ridgepack_timeconvert(nc.time.data(:),...
                            'days since 0000-01-01 00:00:0.0',...
                            nc.time.calendar,1);

         % take difference to get weightings
         diffs=diff(time);

         % if any diffs greater than a month, stop running
         % if any(diffs>28)
         % error('time seperation is more than 28 days apart')
         if any(diffs>31)
          if all(diffs>31)
           error('Time seperation is more than 31 days apart in all cases')
          elseif diffs(diffs>31)>2*std(diffs) & all(diff(diffs(diffs<=31))==0)
           % this takes into account case where timeseries is, say, daily data
           % for an individual month spread over several years, so that 
           % the exceptional diffs are jumps from one year to the next, while
           % the non-exceptional time diffs are between days in the month(s)
           disp('WARNING: Reconstructing time weightings statistically')
           diffs(diffs>31)=mean(diffs(diffs(:)<31));
          elseif diffs(diffs>31)>2*std(diffs) & any(diffs==0)
           disp('WARNING: Reconstructing time weightings, some times are identical')
           diffs(diffs>31 & diffs==0)=mean(diffs(diffs(:)<31 & diffs(:)==0));
          else
           error('Unable to reconstruct time weightings')
          end
         else
          disp('Creating linear time weighting with boundary condition')
         end
         
         % calculate weights in the interior of the dataset based on 
         % half the difference between the sample before the current sample,
         % and half after the current sample.
         weight(2:l-1)=(diffs(1:end-1)+diffs(2:end))./2;

         % calculate weights at the boundaries (first and last sample)
         weight(1)=diffs(1);
         weight(l)=diffs(end);

        else

         disp(['Creating linear weightings for ',name])

	 % generate split differences linear time weightings
         diffs=diff(nc.(name).data);

         % calculate weights in the interior
         weight(2:l-1)=(diffs(1:end-1)+diffs(2:end))./2;

         % calculate weights at the boundaries
         weight(1)=diffs(1);
         weight(l)=diffs(end);

        end

       else % if l==1

        weight=1;

       end

       % normalize the weights
       sumweight=sum(weight(:));
       weight(1:l)=weight(1:l)./sumweight;

       % assign to structure
       nc.(name).weight=reshape(weight,size(nc.(name).data));

       if debug; disp(['Finished creating weightings for ',name]); end

   end

  end

 end

 % Apply linear weightings to var dimensions not included as coordinate dimensions.
 % In other words, apply the weightings calculated in the previous loop for each 
 % dimension over which a mean is being taken to each variable that has that
 % dimension.
 for k=1:length(dims)

  dimension=char(dims{k});
  nc=ridgepack_shuffle(nc,{dimension});

  if any(strcmpi(dimension,dimnames))

   for i=1:length(vars)

    name=char(vars{i});

     if any(strcmp(dimension,nc.(name).dimension))

	disp(['Adding weighting for ',dimension,' to ',name]);

	snc=size(nc.(name).data);

	if size(nc.(name).weight)~=snc
	 disp(['  data array size: ',num2str(snc)]);
	 disp(['weight array size: ',num2str(size(nc.(name).weight))]);
	 error(['2: weight and data arrays different size for ',name]);
	end
	newsize=[prod(snc(1:end-1)),snc(end)];
	nc.(name).weight=reshape(nc.(name).weight,newsize);

	if length(nc.(name).dimension)==1
	 for j=1:length(nc.(dimension).weight(:))
	  nc.(name).weight(j)=nc.(name).weight(j).*nc.(dimension).weight(j);
	 end
        else
	 for j=1:length(nc.(dimension).weight(:))
	  nc.(name).weight(:,j)=nc.(name).weight(:,j).*nc.(dimension).weight(j)';
	 end
        end

	nc.(name).weight=reshape(nc.(name).weight,snc);

        nc.(name).weight(isnan(nc.(name).data))=0;

	if debug; disp(['Finished adding weighting for ',dimension,' to ',name]); end

     end

   end

  end

 end

end % end non-fasttrack items

% apply mask to weightings if a mask exists in the nc structure
% by filling weight with NaNs where mask is zero.

masked=cell(length(vars),1);
if isfield(nc,'mask')

  % Replacing mask NaNs
  if any(isnan(nc.mask.data))
   disp('Setting NaN mask values to zero')
   nc.mask.data(isnan(nc.mask.data))=0;
  end

  % check mask is only 1s and 0s
  if any(nc.mask.data~=0 & nc.mask.data~=1)
   error('mask contains values that are neither unity nor zero')
  end

  nc=ridgepack_shuffle(nc,nc.mask.dimension);
  nc.mask.data(nc.mask.data==0)=NaN;

  for i=1:length(vars)
	name=char(vars{i});

	if ~isempty(intersect(nc.(name).dimension,nc.mask.dimension))
 	  if debug; disp(['Applying mask to ',name]); end
	  if isempty(setdiff(nc.(name).dimension,nc.mask.dimension))
           nc.(name).weight=nc.(name).weight.*nc.mask.data;
	  else
	   masksize=size(nc.mask.data);
	   nmdims=ndims(nc.mask.data);
	   datsize=size(nc.(name).data);
	   if length(nc.(name).dimension)>length(datsize)
	    datsize=[datsize ones(1,length(nc.(name).dimension)-length(datsize))];
	   end
	   firstdims=prod(datsize(1:end-nmdims));
	   lastdims=prod(masksize);
	   nc.(name).weight=reshape(nc.(name).weight,[firstdims,lastdims]);
	   for j=1:firstdims
	    nc.(name).weight(j,:)=nc.(name).weight(j,:).*nc.mask.data(:)';
	   end
	   nc.(name).weight=reshape(nc.(name).weight,datsize);
	  end
	  if ~fasttrack; masked{i}={' (masked)'}; end
	  if debug; disp(['Finished applying mask to ',name]); end
	else
	  if debug; disp(['Mask not applicable to ',name]); end
	end
  end
  nc.mask.data(isnan(nc.mask.data))=0;

else
  if debug; disp('This netcdf structure does not contain a mask field'); end
end


% Permute dimensions in nc so that dimensions to be removed are last
if debug; disp('Shuffling dimensions'); tic; end; 
nc=ridgepack_shuffle(nc,dimnames);
if debug; disp('Finished shuffling dimensions'); toc; end; 

% now pass through each variable to be reduced
for i=1:length(vars)

	name=char(vars{i});

	% check to see that this variable will not be removed completely
	if isempty(setdiff(nc.(name).dimension,dimnames))
         disp([name,' removed from the reduced structure']);
	 continue;
        end

	% Get the reduction template
	if debug; disp('Getting the reduction template'); tic; end; 
	[dimremove,template,tempsize,endsize,ant,time]=...
                       ridgepack_varbounds(nc,dimnames,bounds,name);
	if debug; disp('Got the reduction template'); toc; end; 

	if length(dimremove)==length(nc.(name).dimension)

	 % remove variable completely if needed
	 if debug; disp(['Removing the variable ',name]); end
       
        else

	 % Check to see if a slice or mean is being taken based on the number of 
	 % non-NaN elements in the template.
	 if isempty(template(~isnan(template)))

          if debug; disp('This variable will not be reduced'); end

	  ncr.(name)=nc.(name);
	  continue;

	 % take a slice or provide the cumulative mask
	 elseif strcmp('mask',name) | length(find(~isnan(template(:))))==1 ; 
	  varmean=false;
	  meanonly=false;
	  
         % calculate the standard deviation and potentially output sample field
	 else
	  varstd=[name,'_std'];
	  meanonly=true;

         end

	 % Calculate weighting array only valid for dimensions being reduced.
	 weight=squeeze(reshape(nc.(name).weight,[prod(endsize),prod(tempsize)]));

         % The following code in effect does what is provided here, but more efficiently
	 % for j=1:prod(endsize)
	 %  weight(j,:)=squeeze(reshape(template,[1 prod(tempsize)]).*weight(j,:));
         % end
	 xxx=squeeze(reshape(template,[1 prod(tempsize)]));
         for j=1:length(xxx)
	  weight(:,j)=xxx(j)*weight(:,j);
         end

	 % Remove weight arrays from the netcdf stucture now they have been used.
	 nc.(name)=rmfield(nc.(name),'weight');

	 % Reshape the arrays for calculations into 2D arrays, the first dimension
	 % representing all initial dimensions that are to remain, the second dimension
	 % representing all initial dimensions to be averaged. 
         % Create arrays for the weighted mean and weighted standard deviations.
	 nc.(name).data=squeeze(reshape(nc.(name).data,[prod(endsize),prod(tempsize)]));

	 % Calculate weighted mean or slice and write to new netcdf file
	 ncr.(name).data=ones([1 prod(endsize)]); 
	 ncr.(name).data(:)=nansum(nc.(name).data(:,:).*...
                            weight(:,:),2)./(nansum(weight(:,:),2));

	 if varmean;

	  % Calculate weighted standard deviation and sample size
	  disp(['Calculating weighted standard deviation for ',name]);
	  ncr.(varstd).data=zeros([1 prod(endsize)]); 
	  coeff=zeros(endsize); 

	  for j=1:prod(endsize)
           coeff(j)=sum(abs(weight(j,:))>0);
           if coeff(j)>1 & any(abs(nc.(name).data(j,:))>0)
            ncr.(varstd).data(j)=nansum(weight(j,:).*(nc.(name).data(j,:) - ...
                                ncr.(name).data(j)).^2,2);
            ncr.(varstd).data(j)=ncr.(varstd).data(j)/...
                                ((coeff(j)-1)*nansum(weight(j,:),2)./coeff(j));
            ncr.(varstd).data(j)=sqrt(ncr.(varstd).data(j));
           end
          end

	  % write number of samples to netcdf structure
	  varsamp=[name,'_samp'];
	  ncr.(varsamp).data=coeff;
	  ncr.(varsamp).data=squeeze(reshape(ncr.(varsamp).data,endsize));
	  ncr.(varsamp).dimension=ridgepack_setdiff(nc.(name).dimension,dimremove);
	  ncr.(varsamp).type='NC_SHORT';
	  ncr.(varsamp).long_name=[nc.(name).long_name,...
                                    ' SAMPLES',ridgepack_cellcat(ant),char(masked{i})];

	  % write standard deviation to netcdf structure
	  ncr.(varstd).long_name=[nc.(name).long_name,...
                                 ' WEIGHTED STDEV',ridgepack_cellcat(ant),char(masked{i})];
	  ncr.(varstd).data=squeeze(reshape(ncr.(varstd).data,endsize));
	  ncr.(varstd).data(ncr.(varsamp).data==1)=NaN; 
	  ncr.(varstd).dimension=ridgepack_setdiff(nc.(name).dimension,dimremove);
	  ncr.(varstd).type='NC_FLOAT';
	  if isfield(nc.(name),'units'); ncr.(varstd).units=nc.(name).units; end

          % Calculate equivalent sample size based on Wilks (2006)
          % "Statistical Methods in the Atmospheric Sciences" with further 
          % information from the National Institute of Standards. Note
          % that when all weights are equal but not 1, the lag-1 autocorrelation
          % coefficient is identical to the off-diagonal terms of
          % r=corrcoef(nc.(name).data(j,1:end-1),nc.(name).data(j,2:end)); 
          % which is a verification of the code is correctly approximating
          % unbiased weighted variance and covariance for lag-1.
         
	  varequiv=[name,'_equiv'];
	  ncr.(varequiv).data=zeros([1 prod(endsize)]); 

	  disp(['Estimating equivalent sample size for ',name]);
	  for j=1:prod(endsize)
            if coeff(j)>1 & any(abs(nc.(name).data(j,:))>0)

             % calculate lag-1 weights
             w1=squeeze(weight(j,1:end-1));
             w2=squeeze(weight(j,2:end));

             % calculate lag-1 timeseries (offset by 1)
             x1=squeeze(nc.(name).data(j,1:end-1));
             x2=squeeze(nc.(name).data(j,2:end));

             % get lag-1 length
             N=length(x1);

             % calculate reduced means             
             mean1=nansum(w1.*x1)./nansum(w1);
             mean2=nansum(w2.*x2)./nansum(w2);

             % calculated weighted covariance
             covar1=w1.*(x1-mean1);
             covar2=w2.*(x2-mean2);
             covariance=nansum(covar1.*covar2)./((N-1)*nansum(w1.*w2)./N);

             % calculate variance 
             var1=nansum(w1.*(x1-mean1).^2)./((N-1)*nansum(w1)./N);
             var2=nansum(w2.*(x2-mean2).^2)./((N-1)*nansum(w2)./N);
              
             % lag-1 autocorrelation coefficient
             r1=covariance/sqrt(var1.*var2);
             
             ncr.(varequiv).data(j)=max(0,coeff(j).*((1-r1)./(1+r1))); 

            else

             ncr.(varequiv).data(j)=NaN;

            end
          end

          % write equivalent sample size 
	  ncr.(varequiv).long_name=[nc.(name).long_name,' WEIGHTED EQUIVALENT SAMPLE SIZE',...
                                    ridgepack_cellcat(ant),char(masked{i})];
	  ncr.(varequiv).data=squeeze(reshape(ncr.(varequiv).data,endsize));
	  ncr.(varequiv).data(ncr.(varsamp).data==1)=NaN; 
	  ncr.(varequiv).dimension=ridgepack_setdiff(nc.(name).dimension,dimremove);
	  ncr.(varequiv).type='NC_FLOAT';

         end

	 % complete writing the slice or mean output
	 if varmean | meanonly
	   ncr.(name).long_name=[nc.(name).long_name,' WEIGHTED MEAN',...
                                 ridgepack_cellcat(ant),char(masked{i})];
	 else
	   ncr.(name).long_name=[nc.(name).long_name,' ',ridgepack_cellcat(ant),char(masked{i})];
         end
	 ncr.(name).data=squeeze(reshape(ncr.(name).data,endsize));
	 ncr.(name).dimension=ridgepack_setdiff(nc.(name).dimension,dimremove);
	 ncr.(name).type='NC_FLOAT';
	 if isfield(nc.(name),'units'); 
            ncr.(name).units=nc.(name).units; 
         end
 
 	 % reshape the initial data array back to its initial size
 	 nc.(name).data=squeeze(reshape(nc.(name).data,[endsize,tempsize]));

        end
end

% check that any variables exist
if exist('ncr')==0
 error('All variables will be removed using this dimension specifications') 
end

% now reduce coordinate variables

if ~isempty(setdiff(coor,dimc));

 for i=1:length(coor)
  [dimremove,template,tempsize,endsize,ant]=ridgepack_varbounds(nc,dimnames,bounds,char(coor{i}));
  if length(dimremove)==0; unchange(i)=true; else; unchange(i)=false; end
 end

 if any(unchange)
	 for i=1:length(coor)
	  if unchange(i)
	   name=char(coor{i});
	   ncr.(name)=nc.(name);
          end
	 end
 else

	 % find the latitude and longitude variables in the dataset
	 [lonvar,latvar]=ridgepack_findlatlon(nc,coor);

	 if ~isempty(setdiff(coor,{lonvar,latvar}))
	  error(['Don''t know how to reduce coordinate variables: ',...
	         ridgepack_setdiff(coor,{lonvar,latvar})])
	 end

	 % remove them completely if required
         if length(nc.(lonvar).dimension)==length(dimremove)

	  disp('Coordinate variables removed')

	  if isfield(ncr,lonvar); ncr=rmfield(ncr,lonvar); end
	  if isfield(ncr,latvar); ncr=rmfield(ncr,latvar); end

         else

	  template=reshape(template,[1 prod(tempsize)]);

	  % assign coordinate variables
	  co.(lonvar).data=squeeze(reshape(nc.(lonvar).data,[prod(endsize),prod(tempsize)]));
	  co.(latvar).data=squeeze(reshape(nc.(latvar).data,[prod(endsize),prod(tempsize)]));

	  % obtain lat/lon data only where the template allows
	  co.(lonvar).data=co.(lonvar).data(:,~isnan(template));
	  co.(latvar).data=co.(latvar).data(:,~isnan(template));

	  % generate new coordinates using the meanm command

	  for j=1:prod(endsize)
	   [ncr.(lonvar).data(j),ncr.(latvar).data(j)]=...
	    meanm(co.(lonvar).data(j,:),co.(latvar).data(j,:));
          end

	  % reshape final arrays and write to the output netcdf file
	  ncr.(lonvar).data=reshape(ncr.(lonvar).data,endsize);
	  ncr.(lonvar).dimension=ridgepack_setdiff(nc.(lonvar).dimension,dimremove);
	  ncr.(lonvar).long_name=['Mean longitude along track ',ridgepack_cellcat(ant)];
	  ncr.(lonvar).type='NC_FLOAT';
	  if isfield(nc.(lonvar),'units'); ncr.(lonvar).units=nc.(lonvar).units; end

	  ncr.(latvar).data=reshape(ncr.(latvar).data,endsize);
	  ncr.(latvar).dimension=ridgepack_setdiff(nc.(latvar).dimension,dimremove);
	  ncr.(latvar).long_name=['Mean latitude along track ',ridgepack_cellcat(ant)];
	  ncr.(latvar).type='NC_FLOAT';
	  if isfield(nc.(latvar),'units'); ncr.(latvar).units=nc.(latvar).units; end

         end

 end

end

% keep supporting variables that have had not had their dimensions changed
for i=1:length(supp)
 supname=char(supp{i});
 if isempty(intersect(nc.(supname).dimension,dimnames))
  ncr.(supname)=nc.(supname);
  if debug; disp(['Appending supporting variable ',supname,' to the new structure.']); end
 else
  if debug; disp(['Discarding ',supname,' from the new structure.']); end
 end
end

% add dimension variables to ncr for dimensions to removed during reduction
ncrdims=ridgepack_setdiff(ridgepack_union(dimc,dims),dimnames);
for i=1:length(ncrdims)
 if debug; disp(['Adding ',char(ncrdims{i}),' as a dimension to the new nc structure']); end
 ncr.(char(ncrdims{i}))=nc.(char(ncrdims{i}));
end

% complete writing the new netcdf structure
ncr.attributes=nc.attributes;
comment=['Reduced dataset via slices or means through:',ridgepack_cellcat(dimnames)];
if isfield(ncr.attributes,'comment')
 ncr.attributes.comment=[ncr.attributes.comment,'\n ',comment];
else
 ncr.attributes.comment=comment;
end
if isempty(strfind(ncr.attributes.title,' reduced dataset'))
 ncr.attributes.title=[ncr.attributes.title,' reduced dataset'];
end

% groom new structure
[ncr,out]=ridgepack_struct(ncr);

if debug; disp(['...Leaving ',mfilename]); end

