function [nc]=ridgepack_pcolorm(nc,Z,dimnames,bounds,cont,loglin,ref,horiz,colors,colvals)

% ridgepack_pcolorm - Color fills data from an nc structure on map axes for curvilinear coordinates
%
% function [nc]=ridgepack_pcolorm(nc,Z,dimnames,bounds,cont,loglin,ref,horiz,colors)
%
% This function provides styled pcolor plots on a map setup using the
% matlab mapping toolbox. A graded colorbar is provided.  The data must be in curvilinear
% or rectilinear coordinates.  Geodesic meshes are not allowed for this routine.
%
% Inputs:
% nc     - a netcdf structure (see ridgepack_struct for more details). The
%          structure must contain the fields "latitude" and "longitude". 
%
% Z      - character variable naming the data to be shaded. Several
%          calculation options are possible for Z, including _std for 
%          the standard deviation over the bounds provide, or _samp for 
%          the number of samples used in these calculations. See
%          ridgepack_reduce or ridgepack_select for more information.  If Z is set
%          to 'mask' and a mask field exists in the nc structure, the
%          function shades the mask light grey instead of color rendering
%          a field.  No colorbar is added in this case, but this 
%          can be used to over plot an existing field with a grey mask.
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
% cont   - contour range entered as a vector [C1,C2,...,CX]. This may
%          be entered as an empty vector [] to get ncpcolor to choose
%          the contour interval for you, or it may be omitted if
%          no other proceeding arguments are required.
%
% loglin - 'linear' for linear color scaling, or for log scaling the
%          minumum order of magnitude to be plotted. When using this 
%          log option, you need to make three specifications.  Firstly
%          you need to specify the contour range differently than for
%          for a linear plots.  The numbers in the contour range correspond
%          to order of magnitude multipled by the base order of magnitude
%          which is entered in as loglin.  The sign of the contour values
%          indicates the actual sign of the contour bounds.  For example,
%          say you want to specify a range of -1.e-7 to 1.e-3 on a log10
%          scale in intervals of power of 10, with the smallest value
%          represented being +/- 10-9, then you would enter:
%          cont=[-7:1:3] and loglin=10^-9. 
%          You can adjust the reference value to emphasize color contrast
%          where it is most desired.  If, in the above example, you would
%          like to adjust the spacing to indicate minor logarithmic tick
%          marks on the color scale, you would change cont to:
%          cont=[-7:0.1:3] and loglin=10^-9.
%
% ref    - reference point about which the color is centered, and, when
%          the minimum or maximum is equal to ref, all points below or 
%          abover this value, respectively, are not filled with a color.
%          {zero is the default value}
%
% horiz  - 'horizontal' for colorbar along the bottom, 'vertical' for
%          colorbar on the right side of the plot. Enter 'none' if no
%          colorbar is required. {vertical is default}
%
% colors - This sets the color scheme required:
%	   'bluered'   - highest contrast blue fading through red at ref value
%          'greenred'  - green fading through red at ref value
%          'bluegreen' - blue fading through green at ref value
%          'jet'       - standard matlab 'jet' colorscheme.  Ref has no effect.
%          'cool'      - standard matlab 'cool' color scheme.  Ref has no effect.
%          'gray'      - straight gray color scale. Ref has no effect.
%          'parula'    - standard matlab 'parula' color scheme. Ref has no effect.
%          {bluered is the default}
%
% colvals- Cell array of text for each tick mark on the colorbar if required to 
%          be different from the supplied contour values {optional}. 
%	   See ridgepack_colorbar if further explanation is required.
%
% There are two compulsory inputs: nc and Z. The rest may be 
% included, one at a time, as they are needed in the sequence 
% provided here.
%
%
% OUTPUT:
%
%
% nc - nc structure with information added ont he position of 
%      grid cells' centers if the data is plotted as pixelated grid cells.
%      By outputing this information, it can be re-used in future
%      plotting of the data, drastically reducing the time it takes
%      to produce pixellated data on a map.
%
% Ridgepack Version 1.0
% Andrew Roberts, Naval Postgraduate School, March 2018 (afrobert@nps.edu)
%

global debug;
if debug; disp(['Entering ',mfilename,'...']); end

% get current axes
hmap=get(gcf,'CurrentAxes');

% set defaults and check inputs
if nargin<2; 
 error('specify nc and Z'); 
elseif ~isstruct(nc)
 nc
 error('No netcdf structure input')
elseif ~ischar(Z)
 Z
 error('Z must be character inputs')
elseif ~isfield(nc,Z)
 error([Z,' is missing from this dataset']);
end

if nargin<3; 
 dimnames={}; 
elseif ~iscell(dimnames)
 dimnames
 error('Dimnames must be a cell array')
end

if nargin<4;
 bounds={}; 
elseif ~iscell(bounds)
 bounds
 error('bounds must be a cell array')
end

if nargin<6;
 loglin='linear' ; 
elseif isempty(loglin)
 loglin='';
elseif ~(ischar(loglin) | isnumeric(loglin))
 loglin
 error('loglin must be a character input')
end

if nargin<7;
 ref=0.0; 
elseif ~isnumeric(ref) || length(ref)~=1
 ref
 error('ref must be a single number')
end

if nargin<8;
 horiz='vertical'; 
elseif ~ischar(horiz)
 horiz
 error('horiz must be a character input')
end

if nargin<9;
 colors='bluered' ; 
elseif ~ischar(colors)
 colors
 error('colors must be a character input')
end

% check this function can be used on the chosen map
ht=get(gcf,'CurrentAxes');
if ~ismap(ht)
 error('Current axes must be a map')
end
figout=get(ht,'OuterPosition');
fontsize=min(11,max(8,10*sqrt((figout(3).^2)+(figout(4).^2))/sqrt(2)));
set(ht,'fontsize',fontsize);

% get data slice
[nc]=ridgepack_select(nc,'latitude','longitude',Z,dimnames,bounds);
[nc]=ridgepack_gridsquare(nc);

% change units if need be
[z,nc]=ridgepack_standardunits(nc,Z);
[z,nc]=ridgepack_polarplot(nc,Z);
z=squeeze(z);

% get information from the structure
if isfield(nc,'latitude')
 a=squeeze(nc.latitude.data);
else
 error('latitude is not in the structure');
end

if isfield(nc,'longitude')
 b=squeeze(nc.longitude.data);
else
 error('longitude is not in the structure');
end

% get units
if isfield(nc.(Z),'units')
 units=ridgepack_units(nc,Z);
else
 units=[];
end

zz=z;

% Center pcolor cells on grid cell by splitting the grid.
if debug; disp('Pixelating grid cells'); end
nc=ridgepack_gridsquare(nc);
c=nc.latitude_corner.data;
d=nc.longitude_corner.data;
z=ones(size(zz,1)+1,size(zz,2)+1);
z(1:size(zz,1),1:size(zz,2))=zz(:,:);
z(size(zz,1)+1,1:size(zz,2))=zz(size(zz,1),:);
z(1:size(zz,1),size(zz,2)+1)=zz(:,size(zz,2));
z(size(zz,1)+1,size(zz,2)+1)=zz(size(zz,1),size(zz,2));
if isfield(nc,'mask') && all(size(nc.mask.data)==size(zz))
  mask=ones(size(z));
  mask(2:end,2:end)=nc.mask.data;
  mask(1:end-1,1:end-1)=nc.mask.data;
else
  mask=ones(size(z));
end


% Log or linear color scale minumum, plotted value is minvl
% else get contour values automatically if cont not suppplied 
if isnumeric(loglin)

 if loglin==0; error('loglin must be greater than zero for logarithmic scaling'); end

 if min(diff(cont))<1

  if any(cont>0)
   lcont=unique(round(cont(cont>=0)));
   k=0;
   for i=[min(lcont)-1 lcont]
   for j=1:9;
    k=k+1;
    pcont(k)=loglin*j*(10^(abs(i)-1));
   end
   end
   mincont=loglin*10^(abs(min(cont(cont>=0)))-1);
   if any(cont<0); mincont=floor(mincont); end
   maxcont=loglin*10^(abs(max(cont(cont>=0)))-1);
   pcont=[mincont pcont(pcont>mincont & pcont<maxcont) maxcont];
  else
   pcont=[];
  end

  if any(cont<0)
   lcont=unique(round(cont(cont<=0)));
   k=0;
   for i=[lcont min(abs(lcont))-1]
   for j=9:-1:1;
    k=k+1;
    ncont(k)=-loglin*j*(10^(abs(i)-1));
   end
   end
   mincont=-loglin*10^(abs(min(cont(cont<=0)))-1);
   maxcont=-loglin*10^(abs(max(cont(cont<=0)))-1);
   if any(cont>0); maxcont=ceil(maxcont); end
   ncont=[mincont ncont(ncont>mincont & ncont<maxcont) maxcont];
  else
   ncont=[];
  end

  cont=[ncont pcont];

 else

  cont=unique(round(cont));
  cont=sign(cont).*loglin.*10.^(abs(cont)-1);

 end

 cont=cont(cont~=0); % remove zeros from log scale
 z(abs(z)<min(abs(cont))+eps)=sign(z(abs(z)<min(abs(cont))+eps))*min(abs(cont))+eps;

elseif nargin<5;

 cont=ridgepack_contlev(z);

elseif isempty(cont)

 cont=ridgepack_contlev(z);

elseif length(cont)<2

 disp('Contour levels should be at least two elements long');
 c1=cont(1);
 cont=zeros(1, 3);
 cont(1)=c1-eps;
 cont(2)=c1;
 cont(3)=c1+eps;

end

% blockout values related to the reference min or max
if min(cont)==ref 
	blankout=true;
	z(z<ref)=NaN;
elseif max(cont)==ref
	blankout=true;
	z(z>ref)=NaN;
else
	blankout=false;
	zi=z; 
	zi(z<ref)=NaN;
end

% set colorscheme
if strcmp(loglin,'linear')
 ridgepack_colormap(cont,ref,colors);
else
 ridgepack_colormap(cont,ref,colors,true);
end


% get true color array
[zindex,truecolor]=ridgepack_colorindex(z',cont,ref,mask');


% get map object
mapobj=gcm;


% generate patch arrays
imax=size(c,1)-1;
px=ones(4,imax);
py=ones(4,imax);
pz=ones(4,imax);
pc=ones(1,imax,3);


% for 3D plot, need to get altitude of data (set to zero)
if strmatch(mapobj.mapprojection,'globe') 

 if debug; disp('Plotting on 3D globe'); end

 % set altitude before transforming data
 alt=zeros(size(c));
 [c,d,e] = mfwdtran(c,d,alt,gca);

else

 if debug; disp('Plotting on 2D map projection'); end
 mstruct=ridgepack_extendframe;
 [c,d] = mfwdtran(mstruct,c,d,gca,'surface');

 % set altitude after transforming data to be just above the surface
 e=zeros(size(c));

end

for j=1:size(c,2)-1

 % plot row by row of patches in true color
 px=[c(1:imax,j)';c(2:imax+1,j)';c(2:imax+1,j+1)';c(1:imax,j+1)'];
 py=[d(1:imax,j)';d(2:imax+1,j)';d(2:imax+1,j+1)';d(1:imax,j+1)'];
 pz=[e(1:imax,j)';e(2:imax+1,j)';e(2:imax+1,j+1)';e(1:imax,j+1)'];

 pc=truecolor(j,1:imax,:);

 patch(px,py,pz,pc,'EdgeColor','none')

end

% set up color reference line
if isfield(nc,'mask'); zz(nc.mask.data==0)=NaN; end

% find grid type
gridtype=ridgepack_gridtype(nc);

% zip up global plots along x or latitude for mask plotting
if strcmp(gridtype,'global') 
  b=wrapTo360(b);
  if find(strcmp('x',nc.longitude.dimension))==1 & length(nc.longitude.dimension)>1
   a(end+1,:)=a(1,:);
   b(end+1,:)=360+b(1,:);
   zz(end+1,:)=zz(1,:);
  elseif find(strcmp('x',nc.longitude.dimension))==2 & length(nc.longitude.dimension)>1
   a(:,end+1)=a(:,1);
   b(:,end+1)=360+b(:,1);
   zz(:,end+1)=zz(:,1);
  elseif find(strcmp('longitude',nc.(Z).dimension))==1 
   b(end+1)=360+b(1);
   zz(end+1,:)=zz(1,:);
  elseif find(strcmp('longitude',nc.(Z).dimension))==2 
   b(end+1)=360+b(1);
   zz(:,end+1,:)=zz(:,1);
  end
end

% add dividing line along reference contour unless using log shading
if not(isempty(cont(cont==ref))) & not(blankout) & ~isnumeric(loglin)
 zz(zz<ref)=0;
 zz(zz>=ref)=1;
 ridgepack_maskm(a,b,zz,0.7*[0.99 0.99 0.99]); 
end

% plot mask outline over the top (regardless of if mask is true or not)
%if isfield(nc,'mask')
% zz(~isnan(zz))=1;
% zz(isnan(zz))=0;
% ridgepack_maskm(a,b,zz,[0.5 0.5 0.5]);
%end

% generate colorbar
if strcmp(horiz,'none')
 disp('No colorbar being plotted')
elseif nargin==10 
 ridgepack_colorbar(cont,units,loglin,horiz,ref,colvals); 
else
 ridgepack_colorbar(cont,units,loglin,horiz,ref); 
end

% add title, removing previous title
ridgepack_title(nc,['Shading: ',char(nc.(Z).long_name)],1); 

drawnow;

% clear output if none requested
if nargout==0; clear nc; end

if debug; disp(['...Leaving ',mfilename]); end

