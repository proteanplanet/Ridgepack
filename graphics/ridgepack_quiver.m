function [nc]=ridgepack_quiver(nc1,X,Y,U,nc2,V,dimnames,bounds,thin,mode,scale,underlay,cont)

% ridgepack_quiver - Plot 2D vectors from an nc structure on Cartesian axes
%
% function [nc]=ridgepack_quiver(nc1,X,Y,U,nc2,V,dimnames,bounds,thin,mode,scale,underlay,cont)
%
% This function provides a styled quiver plot on a catesian axis.
%
% INPUT:
%
% nc1    - a netcdf structure for u (see ridgepack_struct for more details)
%
% X      - character variable naming the x-direction coordinate in nc1
%
% Y      - character variable naming the y-direction coordinate in nc1
%
% U      - character variable naming the u-component data 
%
% nc2    - a netcdf structure for v (see ridgepack_struct for more details)
%
% V      - character variable naming the v-component data 
%
% dimnames - cell array of the dimensions to be removed from the nc
%          structure to obtain the requested cross section. See example
%          provided under ridgepack_reduce for more information.
%
% bounds - cell array of the boundaries on each dimension reduction.
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
% mode   - the mode of plotting vectors.  The default is 'b' for 
%          straight length-referenced vectors blue vectors. This can be
%          set to several matlab colors:
%          'b':blue, 'k':black, 'w':white, 'g':green, 'r':red
%          Alternatively, mode can be entered as the text string 'color', 
%          and each vector will be color coded according to its magnitude. 
%
% scale  - if using 'blue' vectors, this is the value given to the 
%          reference vector.  If using 'color', this is the 
%          color division range (minx:x:maxx) of non-zero values. If omitted, 
%          these contour values are calculated automatically or the
%          reference vector is scaled against the median of the input
%          vector field over the domain of the data being plotted.
%
% underlay - optional request to color-underlay scalar fields derived from
%          the vector field. Possible options are:
%	   'curl' for the vorticity of the field
%          'div' for the field divergence
%          Remove this option is no underlay is required.  This does not
%          for the 'color' mode.
%
% cont   - contour interval of the underlay field if required entered as
%          a vector [C1,C2,...,CX]. This is optional.
%
%
% OUTPUT:
% 
% nc     - nc structure with added scalar fields according to the 'underlay'
%          option as detailed above. If no underlay is added, this is not 
%          generated. The nc structure provides all fields in the nc1 and nc2
%          structure for the slice requested for the plot, in addition
%          to the new field. The output and plotted field is at the original
%          resolution of the data, not the thinned vector field.
%
% There are 6 compulsory inputs: nc1, X, Y, U, nc2 and V. The rest may be 
% included, one at a time, as they are needed in the sequence 
% provided here.
%
% Ridgepack Version 1.0
% Andrew Roberts, Naval Postgraduate School, March 2018 (afrobert@nps.edu)


global debug;
if debug; disp(['Entering ',mfilename,'...']); end

ht=get(gcf,'CurrentAxes');

% set defaults
if nargin<6 ; 
 error('specify nc1, X, Y, U, nc2, V') ; 
elseif ~isstruct(nc1) || ~isstruct(nc2)
 nc1
 nc2
 error('nc1 and nc2 must be nc structures')
elseif ~isfield(nc1,X) || ~isfield(nc1,Y) || ~isfield(nc1,U) || ~isfield(nc2,V)
 error('One of X, Y, U or V is missing from the nc structures')
end

if nargin<7; dimnames={}; end
if nargin<8; bounds={}; end
if nargin<9 ; thin=2 ; end
if nargin<10 ; mode='b' ; end
if nargin<12 ; underlay='none' ; end

% get data slice and check for units change
[nc1]=ridgepack_select(nc1,X,Y,U,dimnames,bounds);
[u,nc1]=ridgepack_standardunits(nc1,U);
if strcmp(inputname(1),inputname(5))
 nc2=nc1;
 v=nc2.(V).data;
else
 [nc2]=ridgepack_select(nc2,X,Y,V,dimnames,bounds);
 [v,nc2]=ridgepack_standardunits(nc2,V);
end


% get information from the structure
a=squeeze(nc1.(X).data);
if size(a,1)>1 & size(a,2)>1; error('X must have a single dimension - 1');end
b=squeeze(nc1.(Y).data);
if size(b,1)>1 & size(b,2)>1; error('Y must have a single dimension - 1');end

% check dimension are identical for both nc1 and nc2
if not(isfield(nc2,X)) || not(isfield(nc2,Y))
	error('No fields for X or Y in nc2');
else
	c=squeeze(nc2.(X).data);
	if size(c,1)>1 & size(c,2)>1; error('X must have a single dimension - 2');end
	d=squeeze(nc2.(Y).data);
	if size(d,1)>1 & size(d,2)>1; error('Y must have a single dimension - 2');end

	if not(size(c) == size(a))
		error('x: nc1 and nc2 dimensions are different')
	end

	if not(size(d) == size(b) )
		error('y: nc1 and nc2 dimensions are different')
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

% check that dimensions have not been switched
% and if so do the transpose
if length(a) == size(u,1) & length(b) == size(u,2) & ...
   strcmp(char(nc1.(U).dimension{1}),X) 	;
	try
         u=u';
	 v=v';
	catch
	 disp('unable to transpose - do you need to specify the time required?')
	 error('Unable to transpose')
        end
end

% generate the underlay field and add it to a new nc struture if required
if strcmp(underlay,'none')
 disp('No underlay requested')
elseif strcmp(mode,'color')
 error('An underlay cannot be requested with color coded vectors')
else
 nc=nc1;
 nc.(V)=nc2.(V);
 if strcmp(underlay,'curl')
  nc.curl.data=curl(nc1.(X).data,nc1.(Y).data,nc1.(U).data,nc2.(V).data);
  nc.curl.long_name=['Vorticity for ',nc.(U).long_name];
  nc.curl.units='/s';
  nc.curl.dimension=nc.(U).dimension;
 elseif strcmp(underlay,'div')
  nc.div.data=curl(nc1.(X).data,nc1.(Y).data,nc1.(U).data,nc2.(V).data);
  nc.div.long_name=['Divergence for ',nc.(U).long_name];
  nc.div.units='/s';
  nc.div.dimension=nc.(U).dimension;
 else
  error('Requested underlay not recognized')
 end
 opaque=ones(size(nc.(underlay).data));
 opaque(isnan(nc.(underlay).data))=0.002;

 if nargin<13; 
  cont=ridgepack_contlev(nc.(underlay).data); 
 elseif isempty(cont)
  cont=ridgepack_contlev(nc.(underlay).data); 
 elseif length(cont)<2
  disp('Contour levels should be at least two elements long');
  c1=cont(1);
  cont=zeros(1, 3);
  cont(1)=c1-eps;
  cont(2)=c1;
  cont(3)=c1+eps;
 end

 zindex=ridgepack_colorindex(nc.(underlay).data,cont);
 image(nc.(X).data,nc.(Y).data,zindex,'AlphaData',opaque)
 uunits=ridgepack_units(nc,underlay);
 nccolor(cont,uunits);
 axis xy
 hold on;
 ridgepack_title(nc,['Shading: ',char(nc.(underlay).long_name)],1);
end

% if units of x and y are identical, then make prepare mesh
[x,y]=meshgrid(a,b);

% thin vectors
if thin>0; 
 [x,y,u,v]=ridgepack_quiverthin(x,y,u,v,thin,3); 
elseif thin<0
 [x,y,u,v]=ridgepack_quiverthin(x,y,u,v,abs(thin),2); 
end

% reset a and b
a=squeeze(y(1,:))'; b=squeeze(x(:,1));

% quiver the data.
if strcmpi(mode,'b') | strcmpi(mode,'r') | strcmpi(mode,'k') | ...
   strcmpi(mode,'w') | strcmpi(mode,'g')
  if nargin<11
   % scale vectors against the median
   ridgepack_quiverref(x,y,u,v,units,'median',mode);
  else
   % scale vectors manually and plot the reference vector
   ridgepack_quiverref(x,y,u,v,units,scale,mode);
  end
elseif strcmpi(mode,'color')
  if nargin<11
   scale=ridgepack_contlev(sqrt(u.^2+v.^2));
  end
  ridgepack_quiverref(x,y,u,v,units,'median','col',scale);
else
  error('mode option is incorrect')
end

% the aspect ratio 1:1.
if not(isfield(nc1.(X),'units') & isfield(nc1.(Y),'units')) || ...
   strcmp(char(nc1.(X).units),char(nc1.(Y).units)) 
  axis equal;
  axis tight;
  set(gca,'Layer','bottom');
end

% Add the mask if it exists (NaN values)
if isfield(nc1,'mask')
 zi=nc1.mask.data; zi(nc1.mask.data==0)=NaN;
 if length(nc1.(X).data(:))==size(zi,2)
  ridgepack_mask(nc1.(X).data(:),nc1.(Y).data(:),zi',0.5*[1 1 1]);
 elseif length(nc1.(X).data(:))==size(zi,1)
  ridgepack_mask(nc1.(X).data(:),nc1.(Y).data(:),zi,0.5*[1 1 1]);
 else
  error('Problem with mask');
 end
end

% label the axes
ridgepack_labelxy(nc1,X,Y);

% reposition colorbar if necessary to account for axes labels
if strcmp(mode,'color')
 ridgepack_cbpos(gca,'vertical');
end

% add title
ridgepack_title(nc1,['Vectors: ',char(nc1.(U).long_name)],1);

drawnow

% remove nc structure if not output is requested
if nargout==0 & ~strcmp(underlay,'none') ; clear nc; end

if debug; disp(['...Leaving ',mfilename]); end


