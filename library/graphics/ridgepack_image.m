function [nc]=ridgepack_image(nc,X,Y,Z,dimnames,bounds,cont,loglin,ref,horiz,colors,colvals)

% ridgepack_image - Color fills data from an nc structure on Cartesian axes
%
% function [nc]=ridgepack_image(nc,X,Y,Z,dimnames,bounds,cont,loglin,ref,horiz,colors,colvals)
%
% This function provides a styled filled color plot with a colorbar. This is 
% an alternative to ncpcolor, and uses imagesc in instead. The advantage of using
% this function is that it gives a higher quality plot, with light gray mask
% and transparent areas not plotted as a result of specifying ref.  The disadvantage
% is that it is a more intensive to plot, taking longer.
%
% INPUTS:
%
% nc     - a netcdf structure (see ridgepack_struct for more details)
%
% X      - character variable naming the x-direction coordinate in nc.
%
% Y      - character variable naming the y-direction coordinate in nc
%
% Z      - character variable naming the data to be shaded.  Several
%          calculation options are possible for Z, including _std for 
%          the standard deviation over the bounds provide, or _samp for 
%          the number of samples used in these calculations. See
%          ridgepack_reduce or ridgepack_select for more information.
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
%          be entered as an empty vector [] to get ridgepack_image to choose
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
% Apart from the first six inputs, the rest are identical to 
% input to nccolor.m, the subordinate of this program. There are 
% four compulsory inputs: nc, X, Y and Z. The rest may be 
% included, one at a time, as they are needed in the sequence 
% provided here.
%
%
% OUTPUT:
%
% nc  - nc structure with information added on weigtings, etc
%       for future calls if a movie is being made from 
%       an nc structure. This includes circumstances where
%       a dimension is removed from an nc structure
%
% Ridgepack Version 1.0
% Andrew Roberts, Naval Postgraduate School, March 2018 (afrobert@nps.edu)


global debug;
if debug; disp(['Entering ',mfilename,'...']); end

% set defaults and check inputs
if nargin<4; 
 error('specify nc, X, Y, Z'); 
elseif ~isstruct(nc)
 nc
 error('No netcdf structure input')
elseif ~ischar(X) || ~ischar(Y) || ~ischar(Z)
 X
 Y
 Z
 error('X, Y and Z must be character inputs')
elseif ~isfield(nc,X)
 error([X,' is missing from this dataset']);
elseif ~isfield(nc,Y)
 error([Y,' is missing from this dataset']);
elseif ~isfield(nc,Z)
 error([Z,' is missing from this dataset']);
elseif ~isfield(nc.(X),'dimension') | ~isfield(nc.(Y),'dimension')
 error([X,' or ',Y,' are missing dimensions'])
elseif ~iscell(nc.(X).dimension) | ~iscell(nc.(Y).dimension)
 error([X,' or ',Y,' are dimension field not cell arrays'])
elseif length(nc.(X).dimension)>1
 disp([X,' dimensions are: ',ridgepack_cellcat(nc.(X).dimension)])
 error([X,' must have a single dimension as the axis of a cartesian plot'])
elseif length(nc.(Y).dimension)>1
 disp([Y,' dimensions are: ',ridgepack_cellcat(nc.(Y).dimension)])
 error([Y,' must have a single dimension as the axis of a cartesian plot'])
end

if nargin<5; 
 dimnames={}; 
elseif ~iscell(dimnames)
 dimnames
 error('Fifth Argument: Dimnames must be a cell array')
end

if nargin<6; 
 bounds={}; 
elseif ~iscell(bounds)
 bounds
 error('Sixth Argument: bounds must be a cell array')
end

if nargin<8; 
 loglin='linear' ; 
elseif ~(ischar(loglin) | isnumeric(loglin))
 loglin
 error('Eighth Argument: loglin must be a character or numeric input')
end

if nargin<9; 
 ref=0.0; 
elseif ~isnumeric(ref) || length(ref)~=1
 ref
 error('Ninth Argument: ref must be a single number')
end

if nargin<10; 
 horiz='vertical'; 
elseif ~ischar(horiz)
 horiz
 error('Tenth Argument: horiz must be a character input')
end

if nargin<11; 
 colors='bluered' ; 
elseif ~ischar(colors)
 colors
 error('Eleventh Argument: colors must be a character input')
end


% get data slice
[nc]=ridgepack_select(nc,X,Y,Z,dimnames,bounds);

% shuffle data to get z in correct order
nc=ridgepack_shuffle(nc,{X Y});

% change units if need be
[z,nc]=ridgepack_standardunits(nc,Z);

% check that there isn't a time dimension
if length(size(z))==3 & size(z,1)==1 & strcmp(char(nc.(Z).dimension{1}),'time') ; 
 timed=true; 
elseif length(size(z))==3 & strcmp(char(nc.(Z).dimension{1}),'time') 
 error('Dimension settings are incorrect - too many time frames')
elseif length(size(z))>3 
 error('Dimension settings are incorrect')
else 
 timed=false; 
end
z=squeeze(z);

% get information from the structure
a=squeeze(nc.(X).data);
if size(a,1)>1 & size(a,2)>1; error('X must have a single dimension');end
b=squeeze(nc.(Y).data);
if size(b,1)>1 & size(b,2)>1; error('Y must have a single dimension');end

% get units
if isfield(nc.(Z),'units')
	units=ridgepack_units(nc,Z);
else
	units=[];
end

% log or linear color scale minimum plotted value is minvl
% else get contour values automatically if cont not supplied 
if isnumeric(loglin) && loglin~=0

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

elseif nargin<7 || isempty(cont) ;

 cont=ridgepack_contlev(z);

elseif length(cont)<2 

 disp('Contour levels should be at least two elements long'); 
 c1=cont(1); 
 cont=zeros(1, 3);
 cont(1)=c1-eps;
 cont(2)=c1;
 cont(3)=c1+eps;

elseif isnumeric(loglin)

 error('Problem with loglin specification')

end

% remove data below or above ref, then draw data
blankout=false;
if min(cont)==ref 
 z(z<ref)=NaN;
 blankout=true;
elseif max(cont)==ref
 z(z>ref)=NaN;
 blankout=true;
end

% build the colormap with this information
if strcmp(loglin,'linear')
 ridgepack_colormap(cont,ref,colors);
else
 ridgepack_colormap(cont,ref,colors,true);
end

% Generate mask array
if isfield(nc,'mask') && all(size(nc.mask.data)==size(z))
 mask=nc.mask.data;
else
 mask=ones(size(z));
end

% create axes on which to plot data if they don't already exist
h=get(gcf,'CurrentAxes');
if isempty(h)
 h=axes('xlim',[min(a(:)) max(a(:))],'ylim',[min(b(:)) max(b(:))]);
end

% set fontsize
figout=get(h,'OuterPosition');
fontsize=min(11,max(8,10*sqrt((figout(3).^2)+(figout(4).^2))/sqrt(2)));
set(h,'fontsize',fontsize);

% set the color indices, taking into account switching of axes if needed
if ~isvector(a) | ~isvector(b)
 error('x and/or y axis is a matrix instead of a vector')
elseif length(a)==size(z,1) & length(b)==size(z,2)
 [zindex,truecol]=ridgepack_colorindex(z,cont,ref,mask);
elseif length(a)==size(z,2) & length(b)==size(z,1)
 disp(['length of a is ',num2str(length(a))])
 disp(['length of b is ',num2str(length(b))])
 disp(['size of z is ',num2str(size(z))])
 error('size of axes and data incorrect');
end

% pixelate coordinates
c=[1:length(a)+1]; d=[1:length(b)+1];
c(1)=a(1); c(2:end-1)=(a(1:end-1)+a(2:end))/2; c(end)=a(end);
d(1)=b(1); d(2:end-1)=(b(1:end-1)+b(2:end))/2; d(end)=b(end);
[c,d]=meshgrid(c,d);

% get sizes of x, y, and color
if debug;
 disp(['      Size of a: ',num2str(size(a))])
 disp(['      Size of b: ',num2str(size(b))])
 disp(['Size of truecol: ',num2str(size(truecol))])
end


% stripe the figure with patches
for j=1:size(c,2)-1
 px=[c(1:end-1,j)';c(2:end,j)';c(2:end,j+1)';c(1:end-1,j+1)'];
 py=[d(1:end-1,j)';d(2:end,j)';d(2:end,j+1)';d(1:end-1,j+1)'];
 pc=truecol(j,1:size(c,1)-1,:);
 patch(px,py,pc,'EdgeColor','none')
end


% add white dividing line at reference contour unless using log shading
if not(isempty(cont(cont==ref))) && not(blankout) && ~isnumeric(loglin)
 ridgepack_mask(a,b,z,'w',ref);
 disp('White divider added')
end


% Draw a black line around the mask perimeter 
if isfield(nc,'mask')
 mask(mask<0.5)=NaN;
 ridgepack_mask(a,b,mask);
 disp('Mask outline added')
end

% keep axes on top
set(gca,'Layer','top');
box on;

% label the axes
ridgepack_labelxy(nc,X,Y);

% build the colorbar
if strcmp(horiz,'none')
 disp('No colorbar being plotted')
elseif nargin==12 
 ridgepack_colorbar(cont,units,loglin,horiz,ref,colvals); 
else
 ridgepack_colorbar(cont,units,loglin,horiz,ref); 
end

% add title, removing previous title
if timed
 ridgepack_title(nc,['Shading: ',char(nc.(Z).long_name),' ',datestr(nc.time.data)],1)
else
 ridgepack_title(nc,['Shading: ',char(nc.(Z).long_name)],1)
end

% if units of x and y are identical, then make
% the axes aspect ratio 1:1.
if (isfield(nc.(X),'units') & isfield(nc.(Y),'units') && ...
   strcmp(char(nc.(X).units),char(nc.(Y).units))) | ...
   (~isfield(nc.(X),'units') & ~isfield(nc.(Y),'units'))
  axis tight
  axis equal;
else
  axis normal;
end

hold off;

drawnow;

% clear output if none requested
if nargout==0; clear nc; end

if debug; disp(['...Leaving ',mfilename]); end

