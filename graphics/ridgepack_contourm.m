function [hpos,href,hneg]=ridgepack_contourm(nc,Z,dimnames,bounds,cont,ref,mode,textlist)

% ridgepack_contourM - Contour data on map axes
%
% function [hpos,href,hneg]=ridgepack_contourm(nc,Z,dimnames,bounds,cont,ref,mode,textlist)
%
% This function provides a styled contours on a map.
%
% Inputs:
% nc     - a netcdf structure (see ridgepack_struct for more details)
%
% Z      - character variable naming the data to be contoured
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
% cont   - contour range (minx:x:maxx). To draw a single contour value 
%          enter [conval conval] for cont, where conval is the single 
%          contour value you wish to contour. This may also be entered 
%          as an empty array [] to allow the function to choose the 
%          contour levels.
%
% ref    - reference point about which the contour is dashed. If the
%          contour interval does not include the reference default
%          reference value (zero), then the reference value is not
%          contoured. {zero is the default value}
%
% mode   - The color and mode of contouring.  The default is 'black'. 
%          Other modes are: 'red', 'blue', 'green' and 'white' for colored 
%          contours, or 'blacknolabel' for black, unlabelled contours. 
%
% textlist List of contours to be labeled.  This is an optional vector
%          of numeric values to be labeled on contours.  If this is
%          omitted, then every second contour is labelled. The reference
%          contour is never labeled.
%
%
% OUTPUT:
%
% hpos   - handle of contours greater than ref
% href   - handle of ref contour
% hneg   - handle of contours less than ref
%
% There are 2 compulsory inputs: nc and Z. The rest may be 
% included, one at a time, as they are needed in the sequence 
% provided here. This does not work with a 3D globe projections.
%
% Ridgepack Version 1.0
% Andrew Roberts, Naval Postgraduate School, March 2018 (afrobert@nps.edu)
%

global debug;
if debug; disp(['Entering ',mfilename,'...']); end

labels=true;

% get current map axes
hmap=get(gcf,'CurrentAxes');
if ~ismap(hmap)
 error('Current axes must be a map')
end

% set default line width
lwid=0.25;

% fontsize of contour labels depends on whether using a multiplot or not
r=getappdata(gcf,'MultiplotConfiguration');
if isempty(r)
 fontsize=8;
else
 fontsize=6;
end

disp('WARNING: Contour labels not working with r2014b')

% set defaults
if nargin<2 ; error('specify nc and Z') ; end
if nargin<3; dimnames={}; end
if nargin<4; bounds={}; end
if nargin<6 | isempty(ref) ; ref=0.0 ; end
if nargin<7 | isempty(mode) ; mode='black' ; end

% get data slice
[nc]=ridgepack_select(nc,'latitude','longitude',Z,dimnames,bounds);

% change units to standard metocean units
[z,nc]=ridgepack_standardunits(nc,Z);

% regrid data if on a lat/lon grid to ssmi and polar stereographic
[z,nc]=ridgepack_polarplot(nc,Z);
z=squeeze(z);

% get information from the structure
if isfield(nc,'latitude')
 a=squeeze(nc.latitude.data);
else
 error('latitude is not in the structure');
end
if isfield(nc,'longitude')
 b=wrapTo360(squeeze(nc.longitude.data));
else
 error('longitude is not in the structure');
end


% get units
if isfield(nc.(Z),'units')
 units=ridgepack_units(nc,Z);
else
 units=[];
end

if isfield(nc,'mask') & all(size(nc.mask.data)==size(z))
 z(nc.mask.data==0)=NaN;
end

% zip up global plots along x or latitude
gridtype=ridgepack_gridtype(nc);
if strcmp(gridtype,'global')
  if find(strcmp('x',nc.longitude.dimension))==1 & length(nc.longitude.dimension)>1
   a(end+1,:)=a(1,:);
   b(end+1,:)=360+b(1,:);
   z(end+1,:)=z(1,:);
  elseif find(strcmp('x',nc.longitude.dimension))==2 & length(nc.longitude.dimension)>1
   a(:,end+1)=a(:,1);
   b(:,end+1)=360+b(:,1);
   z(:,end+1)=z(:,1);
  elseif find(strcmp('longitude',nc.(Z).dimension))==1
   b(end+1)=360+b(1);
   z(end+1,:)=z(1,:);
  elseif find(strcmp('longitude',nc.(Z).dimension))==2
   b(end+1)=360+b(1);
   z(:,end+1,:)=z(:,1);
  end
end

% set up mesh for shading the data
c=a;
d=b;
if (size(a,1)==1 | size(a,2)==1) & (size(b,1)==1 | size(b,2)==1)
 [c d z]=ridgepack_mesh(a,b,z);
elseif (size(a,1)==1 | size(a,2)==1) | (size(b,1)==1 | size(b,2)==1)
 error('lat and long variables are dimensioned differently');
end

% get contour values automatically if cont not suppplied 
if nargin<5 | isempty(cont) ;

  cont=ridgepack_contlev(double(z));

  % to be consistent with ridgepack_pcolorm and ridgepack_image, remove top 
  % and bottom contours if not reference value
  if cont(1)~=ref; cont=cont(2:end); end
  if cont(end)~=ref; cont=cont(1:end-1); end
  
end

% sort contour levels into negative, positive and remove reference line
cont1=[];
cont2=[];
j=0;
k=0;
for i=1:length(cont)
	if cont(i)<ref
		k=k+1;
		cont1(k)=cont(i);
	elseif cont(i)>ref
		j=j+1;
		cont2(j)=cont(i);
	end
end

% check that reference contour should be included at all
if (min(cont(:))>ref & ref==0) | (max(cont(:))<ref & ref==0) | isempty(ref)
 noref=true;
else
 noref=false;
end

% set up axes
hm=gcm;
if strmatch(hm.mapprojection,'globe')
 error('ridgepack_contourm does not work with globe projections')
else
 % get coordinates in mapping coordinates and make allowance for stupid
 % matlab trimming rules which clip just inside instead of just outside
 % of mapping frame.
 h=get(gcf,'CurrentAxes');
 mstruct=ridgepack_extendframe;
 [x,y,alt]=mfwdtran(mstruct,c,d,h,'surface'); % trimmed values
end

% arrange colors
if strcmpi(mode,'black') || strcmpi(mode,'blacknolabel');
 disp('Plotting black contours')
 col='k';
elseif strcmpi(mode,'red')
 disp('Plotting red contours')
 col='r';
elseif strcmpi(mode,'blue')
 disp('Plotting blue contours')
 col='b';
elseif strcmpi(mode,'green')
 disp('Plotting green contours')
 col='g';
elseif strcmpi(mode,'white')
 disp('Plotting white contours')
 col='w';
else
 error('Choosing colors')
end

if debug 
 disp(['Size of x is: ',num2str(size(x))])
 disp(['Size of y is: ',num2str(size(y))])
 disp(['Size of z is: ',num2str(size(z))])
end

% plot contours greater than ref
if ~isempty(cont2)
 [C,hpos]=contour(x,y,z,cont2,'LineColor',col,'LineWidth',lwid);
 if ~strcmp(mode,'blacknolabel')
  if nargin<8
   TextList=get(hpos,'TextList');
   set(hpos,'TextList',TextList(1:2:end))
  elseif isnumeric(textlist) & isvector(textlist)
   set(hpos,'TextList',textlist)
  else
   error('textlist must be a vector of numeric contour labels')
  end
  if labels
   clabel(C,hpos,'LabelSpacing',500,'FontSize',fontsize,'Color',col);
   %clabel(C,'space',500,'FontSize',fontsize,'Color',col);
   %clabel(C,'FontSize',fontsize,'Color',col);
  end
 end
else
 hpos=[];
end

% plot contours less than ref
if ~isempty(cont1)
 [C,hneg]=contour(x,y,z,cont1,'LineColor',col,'LineStyle','--','LineWidth',lwid);
 if ~strcmp(mode,'blacknolabel')
  if nargin<8
   TextList=get(hneg,'TextList');
   set(hneg,'TextList',TextList(end:-2:1)')
  elseif isnumeric(textlist) & isvector(textlist)
   set(hneg,'TextList',textlist)
  else
   error('textlist must be a vector of numeric contour labels')
  end
  if labels
   clabel(C,hneg,'LabelSpacing',500,'FontSize',fontsize,'Color',col);
   %clabel(C,'Space',500,'FontSize',fontsize,'Color',col);
   %clabel(C,'FontSize',fontsize,'Color',col);
  end
 end
else
 hneg=[];
end

% contour along reference
if ~noref
 [C,href]=contour(x,y,z,[ref ref],'LineColor',col,'LineStyle','-','LineWidth',lwid*3);
end

% plot regional domain outline 
%if strcmp(gridtype,'regional')
% gcm;
% disp('Plotting regional grid outline')
% cc=[c(1,:) NaN c(end,:) NaN c(:,1)' NaN c(:,end)'];
% dd=[d(1,:) NaN d(end,:) NaN d(:,1)' NaN d(:,end)'];
% plotm(cc,dd,'Color',0.75*[1 1 1]);
%end

% add title
if noref
 ridgepack_title(nc,['Contour: ',char(nc.(Z).long_name),...
             ' (',units,')'],1)
else
 ridgepack_title(nc,['Contour: ',char(nc.(Z).long_name),...
             ' (reference = ',num2str(ref),' ',units,')'],1)
end

% set current axes
set(gcf,'CurrentAxes',hmap);

% return contour handle for legends if requested
if nargout==0
 clear hpos hneg href
end

drawnow

if debug; disp(['...Leaving ',mfilename]); end

