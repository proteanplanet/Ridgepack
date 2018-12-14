function [ncnew]=ridgepack_regrid(nc,u,v,opmode,rec1,rec2,n,method)

% ridgepack_regrid - Regrids data in nc structures to polar preset lat-long grids
%
% function [ncnew]=ridgepack_regrid(nc,u,v,opmode,rec1,rec2,n)
%
% This function regrids data on lat-long grids to preset limited area grids
% over the north and south poles or grids from limited area opmodels
% in polar regions. It only operates on 2D, time dependent fields. In other
% words, u (and v) must have two spatial dimensions, and optionally a time
% dimension.
%
% INPUT:
%
% nc        - netcdf structure (see ridgepack_struct for more information)
%
% u         - u-component if vector is being interpolated, else z-component
%             This is a character variable.
%
% v         - v-component if vector is being interpolated
%             This is a character variable.
%             If there is no v component, enter this as ''.
%
% opmode      - type of grid to be used. Options include:
%        	opmode=1: Hibler/Roberts ice-tide grid
%        	opmode=2: SSM/I grid (north pole, 50km resolution) 
%        	opmode=3: SSM/I grid (south pole, 50km resolution)
%        	opmode=4: SSM/I grid (north pole, 5km resolution)
%        	opmode=5: SSM/I grid variab resolution (see res below)
%        	opmode=6: RASM-WRF grid
%        	opmode=7: RASM-POP/CICE grid
%             (see ridgepack_gridgen for further options)
%
%	      Alternatively, if an nc structure is supplied
%             in opmode, then the grid in the nc structure is
%             used instead of the above preset grids, however
%             the nc structure supplied must have been generated
%             using ridgepack_gridgen, or include the gridded fields x, y, 
%             latitude, longitude and a mask as is ouput from ridgepack_gridgen.
%
% rec1,rec2 - start and end records in time record to be interpolated. 
%             Can be omitted or left empty for the entire dataset to 
%             be interpolated. This can be entered either as a 
%             datestring or index from the netcdf file or nc struture
%             in the time vector.
%
% n         - number of points to append to the boundaries to the grid
%             for processing through filtering algorithms. This can 
%             be omitted or left empty as [] if none are to be added to
%             the boundaries.  If a filter order is being used, then 
%             this value should be half the filter order.
% 
% method    - method by which interpolation is done. Default is 'linear',
%             see help page on "griddata" for other options.
%
%
% OUTPUT:
%
% ncnew     - regridded netcdf structure
%
% Ridgepack Version 1.0
% Andrew Roberts, Naval Postgraduate School, March 2018 (afrobert@nps.edu)
%

global debug;
if debug; disp(['Entering ',mfilename,'...']); end

% set defaults
if nargin<5; rec1=1; rec2=1; end
if nargin<7; n=[]; end
if nargin<8; method='natural' ; end

% run basic checks on the structure
if ~isstruct(nc); error('nc is not a structure'); end
if ~isfield(nc,'latitude') | ~isfield(nc,'longitude')
 error('latitude or longitude not in dataset');
end

% check that u and v are in nc and pass important tests
if ~ischar(u) | ~ischar(v)
 error('u and/or v are not character variables')
end

if ~isfield(nc,u)
 error(['Unable to find ',u,' in the nc structure to be interpolated'])
elseif ~strcmp(v,'') & ~isfield(nc,v)
 error(['Unable to find ',v,' in the nc structure to be interpolated'])
end

if ~isfield(nc,'latitude') | ~isfield(nc,'latitude')
 error('Missing latitude and/or longitude from the input dataset')
else
 nc.longitude.data=wrapTo180(nc.longitude.data);
end

if length(nc.(u).dimension)>3 & ndims(squeeze(nc.(u).data))<=3
 try 
  nc=ridgepack_reduce(nc,setdiff(nc.(u).dimension,{'time'  'x'  'y'}));
 catch
  disp('ERROR: Please reduce your dataset to time and two spatial dimensions')
  error('Unable to reduce dataset to three dimensions')
 end
elseif ndims(nc.(u).data)>3
 error([u,' has more than three dimensions'])
end

if ~strcmp(v,'') & ndims(nc.(v).data)>3
 error([v,' has more than three dimensions'])
end

% mesh the lats and longs
if size(nc.latitude.data,1)>1 & size(nc.longitude.data,1)>1 & ...
   size(nc.latitude.data,2)>1 & size(nc.longitude.data,2)>1

   disp('Lats and longs are uniquely represented on the mesh')

   latlon1d=false;
   disp('Using matlab function griddata for interpolation');

   % allow for longitude wrapping
   llat=[nc.latitude.data nc.latitude.data nc.latitude.data];
   llon=[nc.longitude.data-360. nc.longitude.data nc.longitude.data+360.];

   % allow for polar singularity for preset grids
   cutlatt=80; % cutoff for adding in polar value
   maxlatt=90; % 90 degrees north or south 
   if isfield(nc,'stereographic_mesh') | isfield(nc,'icepack_mesh')
    if max(nc.latitude.data(:))>cutlatt
      llat=llat(:); llat(end+1:end+4)=[maxlatt maxlatt maxlatt maxlatt];
    elseif min(nc.latitude.data(:))<-cutlatt
      llat=llat(:); llat(end+1:end+4)=-1*[maxlatt maxlatt maxlatt maxlatt];
    else
      error('no values found close to the pole')
    end
    llon=llon(:); llon(end+1:end+4)=[0.0 90.0 180.0 270.0];
   end

else

   disp('Meshing lat and longs from vectors')

   lat=nc.latitude.data;
   lon=nc.longitude.data;
   llat=nc.latitude.data;
   llon=nc.longitude.data;

   latlon1d=true;
   disp('Using matlab function interp2 for interpolation');

   % allow for longitude wrapping
   llon=[llon-360. llon llon+360.];

   cutlat=87; % cutoff for adding in polar value
   maxlat=90; % 90 degrees north or south 
   if llat(1)<maxlat & llat(1)>cutlat
	   llat(2:length(llat)+1)=llat(1:length(llat));
	   llat(1)=maxlat;
   elseif llat(end)<maxlat & llat(end)>cutlat
	   llat(end+1)=maxlat;
   end
   if llat(1)>-maxlat & llat(1)<-cutlat
	   llat(2:length(llat)+1)=llat(1:length(llat));
	   llat(1)=-maxlat;
   elseif llat(end)>-maxlat & llat(end)<-cutlat
	   llat(end+1)=-maxlat;
   end

   [llat,llon]=meshgrid(llat,llon);

end

% get new lats and longs of grid
if isstruct(opmode)
 ncnew=opmode;
 if ~isfield(ncnew,'latitude') | ~isfield(ncnew,'latitude')
  error('Missing latitude and/or longitude from the opmode grid')
 end
 try
  ncnew=ridgepack_shuffle(ncnew,{'y' 'x'});
 catch
  disp('Cannot find x and y dimensions in provided grid')
 end
else
 ncnew=ridgepack_gridgen([],opmode,n);
end

% wrap to 180
ncnew.longitude.data=wrapTo180(ncnew.longitude.data);

% check for time dimension
timepresent=false;
[variablenames,numbervariables]=ridgepack_name(nc);
for m=1:numbervariables
	name=char(variablenames(m))
	if strcmp(name,'time'); 
	       timepresent=true;
	       if ischar(rec1)
		 rec1=ncnearest(nc.time.data,datenum(rec1));
               end
	       if ischar(rec2)
		 rec2=ncnearest(nc.time.data,datenum(rec2));
               end
	       if nargin<6 | (isempty(rec2) & isempty(rec1)); 
		 rec2=length(nc.time.data); 
		 rec1=1;
	       else
	   	 rec2=min(rec2,length(nc.time.data));
	   	 rec1=min(max(rec1,1),rec2);
	       end
	       disp(['rec1=',num2str(rec1),' rec2=',num2str(rec2)]);
	       ncnew.time=nc.time;
	       ncnew.time.data=nc.time.data(rec1:rec2);
	elseif length(nc.(name).dimension)==1 & ...
	       strcmp(name,char(nc.(name).dimension)) & ...
	       (strcmp(name,'x') | strcmp(name,'y')) 
	       dimx='x';
	       dimy='y';
	elseif length(nc.(name).dimension)==1 & ...
	       strcmp(name,char(nc.(name).dimension)) & ...
	       (strcmp(name,'latitude') | strcmp(name,'longitude'))
	       dimx='latitude';
	       dimy='longitude';
        else
               disp('no assignment')
	end
end

% assign multidimensional data to new netcdf structure
for m=1:numbervariables
	name=char(variablenames(m));
	if timepresent && length(nc.(name).dimension)>2  & ...
           max(strcmp(nc.(name).dimension,dimx))>0 & ...
           max(strcmp(nc.(name).dimension,dimy))>0 & ...
	   (strcmp(name,u) | strcmp(name,v))
	   	if isfield(nc.(name),'units'); ncnew.(name).units=nc.(name).units; end
		ncnew.(name).long_name=nc.(name).long_name;
		ncnew.(name).dimension={'time','y','x'};
	elseif length(nc.(name).dimension)>1 & ...
           max(strcmp(nc.(name).dimension,dimx))>0 & ...
           max(strcmp(nc.(name).dimension,dimy))>0 & ...
	   (strcmp(name,u) | strcmp(name,v))
	   	if isfield(nc.(name),'units'); ncnew.(name).units=nc.(name).units; end
		ncnew.(name).long_name=nc.(name).long_name;
		ncnew.(name).dimension={'y','x'};
	end
end


% set up new structure
if isfield(nc.attributes,'title')
	ncnew.attributes.title=['Interpolated ',nc.attributes.title];
else
	ncnew.attributes.title=['Interpolated data'];
end
ncnew.attributes.source=['bi-cubic interpolation from ',nc.attributes.title,' NetCDF data'];

% shuffle so that time is the last dimension
nc=ridgepack_shuffle(nc,{'time'});


% Remove turning angle of non-vector components are being processed,
% otherwise turn vectors to zonal and meridional if a turning angle
% is supplied. If interpolating to a polar azimuthal projection, rotate
% with longitude, so as to avoid creating a singularity at the pole.
if strcmp(v,'') && isfield(ncnew,'turn')
 ncnew=rmfield(ncnew,'turn');
elseif ~strcmp(v,'') & isfield(nc,'turn')
 disp('Turning vectors to geographic coordinates')
 if ~isstruct(opmode) && opmode>0 & opmode<8 
  exth=deg2rad(nc.longitude.data);
 else
  exth=0.;
 end
 if timepresent & max(strcmp(nc.(u).dimension,'time'))>0
  for i=rec1:rec2
   [th,z]=cart2pol(squeeze(nc.(u).data(:,:,i)),squeeze(nc.(v).data(:,:,i)));
   [nc.(u).data(:,:,i),nc.(v).data(:,:,i)]=pol2cart(th-exth+deg2rad(nc.turn.data),z);
  end
 else
  [th,z]=cart2pol(squeeze(nc.(u).data(:,:)),squeeze(nc.(v).data(:,:)));
  [nc.(u).data(:,:),nc.(v).data(:,:)]=pol2cart(th-exth+deg2rad(nc.turn.data),z);
 end
end

% now interpolate the data
[variablenames,numbervariables]=ridgepack_name(ncnew);
for m=1:numbervariables
	name=char(variablenames(m));
	if not(strcmp(name,'longitude') | strcmp(name,'latitude') | ... 
	   strcmp(name,'turn')) & length(ncnew.(name).dimension)>1 & ...
	   (strcmp(name,u) | strcmp(name,v))

	 if strcmp(name,u) 
	   turnu=name;
	 elseif strcmp(name,v) 
	   turnv=name;
	 end

         if timepresent & max(strcmp(nc.(name).dimension,'time'))>0
	   rec=rec1;
	   sec=rec2;
	 else
	   rec=1;
	   sec=1;
	 end

	 ncnew.(name).data=ones((sec-rec+1),size(ncnew.latitude.data,1),size(ncnew.latitude.data,2));

 	 for i=rec:sec
           

           if timepresent  
            z=squeeze(nc.(name).data(:,:,i));
           else
            z=nc.(name).data;
	   end
	
	   if latlon1d 

            % wrap longitudes to account for -360 to 720 degrees.
            lenlon=length(lon);
            if size(z,2)==lenlon
                  z=[z z z]; 
            else
                  z=[z; z; z]; 
            end

            % fill in poles for a global opmodel
            lenlat=length(lat);
	    if lat(1)<maxlat &  lat(1)>cutlat
	     if size(z,2)==lenlat
	          z(:,2:lenlat+1)=z(:,:);
	          z(:,1)=z(:,2);
		  if i==rec ; disp('Adding NH lat extended data to second z dimension start') ; end
	     else
	          z(2:lenlat+1,:)=z(:,:);
	          z(1,:)=z(2,:);
		  if i==rec ; disp('Adding NH lat extended data to first z dimension start') ; end
	     end
	     lenlat=lenlat+1;
	    elseif lat(end)<maxlat & lat(end)>cutlat
	     if size(z,2)==lenlat
	          z(:,lenlat+1)=z(:,lenlat);
		  if i==rec ; disp('Adding NH lat extended data to second z dimension end') ; end
	     else
	          z(lenlat+1,:)=z(lenlat,:);
		  if i==rec ; disp('Adding NH lat extended data to first z dimension end') ; end
	     end
	     lenlat=lenlat+1;
            end

	    if lat(1)>-maxlat &  lat(1)<-cutlat
	     if size(z,2)==lenlat
	          z(:,2:lenlat+1)=z(:,:);
	          z(:,1)=z(:,2);
		  if i==rec ; disp('Adding SH lat extended data to second z dimension start') ; end
	     else
	          z(2:lenlat+1,:)=z(:,:);
	          z(1,:)=z(2,:);
		  if i==rec ; disp('Adding SH lat extended data to first z dimension start') ; end
	     end
	    elseif lat(end)>-maxlat & lat(end)<-cutlat
	     if size(z,2)==lenlat
	          z(:,lenlat+1)=z(:,lenlat);
		  if i==rec ; disp('Adding SH lat extended data to second z dimension end') ; end
	     else
	          z(lenlat+1,:)=z(lenlat,:);
		  if i==rec ; disp('Adding SH lat extended data to first z dimension end') ; end
	     end
            end

           else
	    % to allow for polar points and longitude wrapping
   	    if isfield(nc,'stereographic_mesh') | isfield(nc,'icepack_mesh')
             if max(nc.latitude.data(:))>cutlatt
   	      nintey=find(abs(nc.latitude.data(:)-maxlatt)==min(abs(nc.latitude.data(:)-maxlatt)));
      	     elseif min(nc.latitude.data(:))<-cutlatt
   	      nintey=find(abs(nc.latitude.data(:)+maxlatt)==min(abs(nc.latitude.data(:)+maxlatt)));
             else
              error('no values found close to the pole')
             end
	     nintey=nanmean(z(nintey));
             z=[z z z]; 
	     z=z(:); z(end+1:end+4)=[nintey nintey nintey nintey];
            else
             z=[z z z]; 
            end
	   end

	   if timepresent 
	    disp(['Interpolating record ',num2str(i),' ',...
                  datestr(nc.time.data(i),0),' for ',name]);
	   end

           if debug
	    disp(['Size of lat: ',num2str(size(llat))])
	    disp(['Size of lon: ',num2str(size(llon))])
	    disp(['Size of z: ',num2str(size(z))])
           end

	   if latlon1d
	    disp('1D lat-lon interpolation')
	    y=interp2(llat,llon,z',ncnew.latitude.data,ncnew.longitude.data,method);
           else
            if any(isnan(ncnew.latitude.data(:))) | any(isnan(ncnew.longitude.data(:)))
             error('New latitude/longitude points have NaNs')
            else
	     disp('2D lat-lon interpolation')

             maskm=(~isnan(llon) & ~isnan(llat) & ~isnan(z));

             llon1d=llon(maskm);
             llat1d=llat(maskm);
             z1d=z(maskm);

	     y=griddata(llat1d,llon1d,z1d,...
                        ncnew.latitude.data,ncnew.longitude.data,method);
	     newmask=griddata(llat1d,llon1d,z1d,...
                        ncnew.latitude.data,ncnew.longitude.data,method);
            end
           end

	   if ~isempty(find(isnan(y))); disp('NaNs found in interpolated field'); end

           if isfield(ncnew,'mask') & all(size(y)==size(squeeze(ncnew.mask.data)))
            y(squeeze(ncnew.mask.data)==0)=NaN;
           end

	   if timepresent
	    ncnew.(name).data(i-rec+1,:,:)=y;
           else
	    ncnew.(name).data=y;
           end

         end

	 ncnew.(name).data=ncnew.(name).data;
	 ncnew.(name).type='NC_FLOAT';

        end
end


% turn vectors to local grid if required, including turning back from 
% longitude if using polar azimuthal grids, which is done to avoid
% a singularity at the pole.
if ~strcmp(v,'') 
  if isfield(ncnew,'turn')
   uvel=ncnew.(turnu).data;
   vvel=ncnew.(turnv).data;
  else
   error('no turning angle provided for new vector components')
  end
  (size(uvel)~=size(vvel))
  if any(size(uvel)~=size(vvel))
   error('u and v are not of equal size')
  end
  disp('Turning the vectors for northern hemisphere');
  if ~isstruct(opmode) && opmode>0 & opmode<8 
   exth=deg2rad(ncnew.longitude.data);
  else
   exth=0.;
  end
  
  if timepresent & max(strcmp(nc.(u).dimension,'time'))>0
   for i=1:size(u,1)
    [th,z]=cart2pol(squeeze(uvel(i,:,:)),squeeze(vvel(i,:,:)));
    [ncnew.(turnu).data(i,:,:),ncnew.(turnv).data(i,:,:)]=...
                pol2cart(th+exth-deg2rad(ncnew.turn.data),z);
   end
  else
   [th,z]=cart2pol(squeeze(uvel(:,:)),squeeze(vvel(:,:)));
   [ncnew.(turnu).data(:,:),ncnew.(turnv).data(:,:)]=...
                pol2cart(th+exth-deg2rad(ncnew.turn.data),z);
  end

  % OLD CODE TO ACCOUNT FOR SOUTHERN HEMISPHERE (PROBABLY WRONG!)
  % if median(ncnew.turn.data)<0
  %  disp('Turning the vectors for southern hemisphere');
  %  [th,z]=cart2pol(-u,-v);
  %  [ncnew.(turnu).data,ncnew.(turnv).data]=pol2cart(-th0deg2rad(ncnew.turn.data),z);
  % end

end

ncnew=ridgepack_struct(ncnew);

if debug; disp(['...Leaving ',mfilename]); end

