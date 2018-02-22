function ridgepack_quiverm(nc1,U,nc2,V,dimnames,bounds,thin,mode,scale)

% ridgepack_quiverm - Plot 2D vectors from an nc structure on map axes.
%
% function ridgepack_quiverm(nc1,U,nc2,V,dimnames,bounds,thin,mode,scale)
%
% This function provides a styled quiver plot on a catesian axis.
%
% Inputs:
% nc1    - a netcdf structure for u (see ridgepack_struct for more details)
%
% U      - character variable naming the u-component data 
%
% nc2    - a netcdf structure for v (see ridgepack_struct for more details)
%
% V      - character variable naming the v-component data 
%
% dimnames - Cell array of the dimensions to be removed from the nc
%          structure to obtain the requested cross section. See example
%          provided under ridgepack_reduce for more information.
%
% bounds - Cell array of the boundaries on each dimension reduction.
%          This includes either a low and high index over which a mean
%          or slice is to be taken, or a low and high value representing
%          the range over which the data is to be extracted.  Time
%          values can be used. See example provided for ridgepack_reduce for
%          more information.
%
% thin   - Thinning to be done as a power of two. thin may be set
%          to 1, to reduce the number of vectors by two, 2 to reduce the
%          number of vectors by four, 3 to reduce the number of vectors
%          by eight, and so on, along each side of a square grid mat. 
%          The standard method is to take the median for all values of the
%          given (2^(thin))^2 grid points in the mat, supplying a median
%          when at least half of the grid points are not NaN values. This
%          is done to accommodate the presence of a land mask. If thin
%          is set to negative, a mean is calculated for each square mat
%          instead of a median. This method only supplies a value if all 
%          values in the mat are not NaNs. Default is 2.
%
% mode   - The mode of plotting vectors.  The default is 'b' for 
%          straight length-referenced vectors blue vectors. This can be
%          set to several matlab colors:
%          'b':blue, 'k':black, 'w':white, 'g':green, 'r':red, 'm':magenta
%          Alternatively, mode can be entered as the text string 'color', 
%          and each vector will be color coded according to its magnitude. 
%
% scale  - If setting monocolor ('b,k,w,g,r,m') vectors, this is the value 
%          given to the reference vector. Setting this value to negative
%          will scale the vectors without plotting the reference vector.
%          This is useful for multi-frame plots.  
%          If using 'color', this is the  contour interval (x) or 
%          range (minx:x:maxx). If omitted, these contour values are 
%          calculated automatically or the reference vector is scaled 
%          against the median of the input vector field over the domain of 
%          the data being plotted. 
%
% There are 4 compulsory inputs: nc1, U, nc2 and V. The rest may be 
% included, one at a time, as they are needed in the sequence 
% provided here. u and v are assumed to be int the +longitude and 
% +latitude directions.
%
% Andrew Roberts, Naval Postgraduate School, March 2018  (afrobert@nps.edu)
%

global debug;
if debug; disp(['Entering ',mfilename,'...']); end

ht=get(gcf,'CurrentAxes');

% check this function can be used on the chosen map
if ismap(ht)
 mapobj=gcm;
 if strmatch(mapobj.mapprojection,'globe')
   error('This function does not work with the globe projections')
 end
else
 error('Current axes must be a map')
end

% set defaults
if nargin<4 ; error('specify nc1, U, nc2, V') ; end
if nargin<5; dimnames={}; end
if nargin<6; bounds={}; end
if nargin<7 ; thin=2 ; end
if nargin<8 ; mode='b' ; end

% get data slice and check for units change
[nc1]=ridgepack_select(nc1,'latitude','longitude',U,dimnames,bounds);
[u,nc1]=ridgepack_standardunits(nc1,U);
if strcmp(inputname(1),inputname(3))
 nc2=nc1;
 v=nc2.(V).data;
else
 [nc2]=ridgepack_select(nc2,'latitude','longitude',V,dimnames,bounds);
 [v,nc2]=ridgepack_standardunits(nc2,V);
end

% regrid data if on a lat/lon grid to ssmi and polar stereographic
[u,nc1]=ridgepack_polarplot(nc1,U);
[v,nc2]=ridgepack_polarplot(nc2,V);

% get information from the structure
if isfield(nc1,'latitude')
 a=squeeze(nc1.latitude.data);
else
 error('latitude is not in the structure');
end
if isfield(nc1,'longitude')
 b=squeeze(nc1.longitude.data);
else
 error('longitude is not in the structure');
end

% check dimension are identical for both nc1 and nc2
if not(isfield(nc2,'latitude')) || not(isfield(nc2,'longitude'))
 error('No fields for X or Y in nc2');
else
 c=squeeze(nc2.latitude.data);
 d=squeeze(nc2.longitude.data);
 if not(size(c) == size(a))
	 error('latitude dimensions are not coherent')
 elseif not(size(d) == size(b))
	 error('longitude dimensions are not coherent')
 end
end

% get units
units=[];
if isfield(nc1.(U),'units') 
	if isfield(nc2.(V),'units')
		if strcmp(char(nc1.(U).units),char(nc2.(V).units))
			units=ridgepack_units(nc1,U);
		else
			error('u and v have different units');
		end
	else
		error('u has units but v has none');
	end
elseif isfield(nc2.(V),'units')
	error('v has units but u has none');
end

% turn vectors if need be
if isfield(nc1,'turn')

	disp('Turning vectors')
	[th,z]=cart2pol(u,v);
	[u,v]=pol2cart(th+deg2rad(nc1.turn.data),z);

        % OLD CODE TO ACCOUNT FOR SOUTHERN HEMISPHERE (PROBABLY WRONG!)
	%if median(nc1.turn.data)<0
	%	[th,z]=cart2pol(-u,-v);
	%	[u,v]=pol2cart(-th-deg2rad(nc1.turn.data),z);
	%end

end

% set up mesh for shading the data
x=a;
y=b;
if (size(a,1)==1 | size(a,2)==1) & (size(b,1)==1 | size(b,2)==1)
	 [x y u]=ridgepack_mesh(a,b,u);
	 [x y v]=ridgepack_mesh(a,b,v);
elseif (size(a,1)==1 | size(a,2)==1) | (size(b,1)==1 | size(b,2)==1)
	 error('lat and long variables are dimensioned differently');
end

% thin vectors
if thin>0; 
 [x,y,u,v]=ridgepack_quiverthin(x,y,u,v,thin,3); 
elseif thin<0
 [x,y,u,v]=ridgepack_quiverthin(x,y,u,v,abs(thin),2); 
end

% quiver the data.
if strcmpi(mode,'b') | strcmpi(mode,'r') | strcmpi(mode,'k') | strcmpi(mode,'w') | strcmpi(mode,'g') | strcmpi(mode,'m')
  if nargin<9
   % scale vectors against the median
   ridgepack_quiverref(x,y,u,v,units,'median',mode);
  else
   % scale vectors manually and plot the reference vector
   ridgepack_quiverref(x,y,u,v,units,scale,mode);
  end
elseif strcmpi(mode,'color')
  if nargin<9
   scale=ridgepack_contlev(sqrt(u.^2+v.^2));
  end
  ridgepack_quiverref(x,y,u,v,units,'median','col',scale);
else
 error('mode option is incorrect')
end

% add title
ridgepack_title(nc1,['Vectors: ',char(nc1.(U).long_name)],1);

drawnow

set(gcf,'CurrentAxes',ht);

if debug; disp(['...Leaving ',mfilename]); end


