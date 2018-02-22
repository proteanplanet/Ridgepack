function [hpos,href,hneg]=ridgepack_contour(nc,X,Y,Z,dimnames,bounds,cont,ref,mode,textlist) 

% ridgepack_contour - Contour data on Cartesian axes
%
% function [hpos,href,hneg]=ridgepack_contour(nc,X,Y,Z,dimnames,bounds,cont,ref,mode,textlist) 
%
% This function provides a styled contours on a Cartesian plot.
%
% Inputs:
% nc     - a netcdf structure (see ridgepack_struct for more details)
%
% X      - character variable naming the x-direction coordinate in nc
%
% Y      - character variable naming the y-direction coordinate in nc
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
%	   Finally, 'phase' plots the cyclical phase of a wave in degrees,
%	   designed for showing tidal phases.
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
% There are 4 compulsory inputs: nc, X, Y and Z. The rest may be 
% included, one at a time, as they are needed in the sequence 
% provided here.
%
% Ridgepack Version 1.0
% Andrew Roberts, Naval Postgraduate School, March 2018 (afrobert@nps.edu)
%

global debug;
if debug; disp(['Entering ',mfilename,'...']); end

% line width
lwid=0.4;

% fontsize of contour labels depends on whether using a multiplot or not
r=getappdata(gcf,'MultiplotConfiguration');
if isempty(r)
 fontsize=6;
else
 fontsize=4;
end

% set defaults
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
end

if nargin<5;
 dimnames={};
elseif ~iscell(dimnames)
 dimnames
 error('Dimnames must be a cell array')
end

if nargin<6;
 bounds={};
elseif ~iscell(bounds)
 bounds
 error('bounds must be a cell array')
end

if nargin<7;
 cont=[];
elseif ~isnumeric(cont)
 cont
 error('cont must be numeric')
end

if nargin<8 || isempty(ref)
 ref=0.0;
elseif ~isnumeric(ref) || length(ref)~=1
 ref
 error('ref must be a single number')
end

if nargin<9;
 disp('Plotting black contours')
 col='k'; 
 mode='black';
elseif ~ischar(mode)
 mode
 error('mode must be a character string')
elseif strcmpi(mode,'black') || strcmpi(mode,'blacknolabel'); 
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
 disp('Plotting green contours')
 col='w'; 
elseif strcmpi(mode,'phase')
 col='k'; 
else	
 disp('Plotting black contours')
 error('')
end


% get data slice
[nc]=ridgepack_select(nc,X,Y,Z,dimnames,bounds);

% change units if need be
[z,nc]=ridgepack_standardunits(nc,Z);
z=squeeze(z);

% set obscissa and ordinate information
a=squeeze(nc.(X).data);
if size(a,1)>1 & size(a,2)>1; error('X must have a single dimension');end
b=squeeze(nc.(Y).data);
if size(b,1)>1 & size(b,2)>1; error('Y must have a single dimension');end


% arrange mask
if isfield(nc,'mask') && all(size(nc.mask.data)==size(z))
 mask=nc.mask.data;
 mask(mask<0.5)=NaN;
else
 mask=NaN;
end


% get units
if isfield(nc.(Z),'units')
	units=ridgepack_units(nc,Z);
else
	units=[];
end


% arrange data if mode='phase' to peak at 180 degrees
if strcmpi(mode,'phase') & strcmpi(units,'degrees')

 disp('Setting up phase')

 % use to arrays to plot phase, one between 0 and 360, the other -180 to 180
 % which corresponse to z and zi respectively.  zi is used for plotting 
 % the zero phase contour only.
 z=wrapTo360(z);
 zi=wrapTo180(z);

 % Fill amphidromic points with NaNs to avoid an unsightly contour tangle.
 % To do this, get 2nd order difference of x- and y-phase components
 % and mask out grid cells where it exceeds 90 degrees 2nd order difference
 % in both directions (this may include the 0-360 degree crossover).
 xdif=abs(diff(z,2,1)); ydif=abs(diff(z,2,2));
 xd=zeros(size(z)); yd=zeros(size(z));
 xd(2:size(xdif,1)+1,:)=xdif; yd(:,2:size(ydif,2)+1)=ydif; 
 z(xd>90 & yd>90)=NaN;

 % mask out areas of the arrays that are not used in plotting to avoid 
 % plotting of contour gradients at the boundary between 0 and 360 and
 % -180 and 180 phase, respectively.
 z(z>350 | z<10)=NaN;
 zi(zi<-90 | zi>90)=NaN;

 % no reference contour is plotted for phase
 noref=true;

elseif strcmpi(mode,'phase') 

 error('Units should be degrees for phase plot')

else

 % get contour values automatically if cont not suppplied 
 if nargin<7 | isempty(cont) ;
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
	if cont(i) < ref
		k=k+1;
		cont1(k)=cont(i);
	elseif cont(i) > ref
		j=j+1;
		cont2(j)=cont(i);
	end
 end

 % check that reference contour should be included at all
 if (min(cont(:))>ref & ref==0) | (max(cont(:))<ref & ref==0) 
	noref=true;
	disp('Removing reference contour')
 else
	noref=false;
	disp('Adding reference contour')
 end

end


% contour the data, noting that transpose must be used for contour
if strcmpi(mode,'phase')

	disp('Plotting phase contours')

	% 30 degree contours
	[C,h]=contour(a,b,z',[30:30:330],'LineColor',col,'Parent',gca); %

	set(h,'ShowText','on','TextStep',get(h,'LevelStep')*2);

	% label 90 degree contours
        set(h,'TextList',[90:90:270]) 

	clabel(C,h,'LabelSpacing',500,'Fontsize',6)

        hold on

	[C,h]=contour(a,b,zi',[0 0],'LineStyle','-','LineColor',col,'Parent',gca);

	set(h,'ShowText','on','TextStep',get(h,'LevelStep')*2);

	clabel(C,h,'LabelSpacing',500,'Fontsize',6)
 
else

	% if black and white, dash the zero contour, and make negative grey
	if ~isempty(cont2)
	 [C,hpos]=contour(a,b,z',cont2,'LineColor',col,'LineWidth',lwid);
	 if ~strcmpi(mode,'blacknolabel'); 
	  if nargin<10
	   TextList=get(hpos,'TextList');
	   set(hpos,'TextList',TextList(1:2:end))
	  elseif isnumeric(textlist) & isvector(textlist)
	   set(hpos,'TextList',textlist)
	  else
	   error('textlist must be a vector of numeric contour labels')
	  end
	  clabel(C,hpos,'LabelSpacing',400,'Fontsize',fontsize,'Color',col,'Margin',1);
         end
	else
	 hpos=[];
	end

	hold on

	if ~isempty(cont1)
	 [C,hneg]=contour(a,b,z',cont1,'LineColor',col,'LineStyle','--','LineWidth',lwid);
	 if ~strcmpi(mode,'blacknolabel'); 
	  if nargin<10
	   TextList=get(hneg,'TextList');
	   set(hneg,'TextList',TextList(end:-2:1)')
	  elseif isnumeric(textlist) & isvector(textlist)
	   set(hneg,'TextList',textlist)
	  else
	   error('textlist must be a vector of numeric contour labels')
	  end
	  clabel(C,hneg,'LabelSpacing',500,'Fontsize',fontsize,'Color',col,'Margin',1);
         end
        else
         hneg=[];
	end

	% contour along reference or for phase contour along 0 (aliased as 180 in zi)
	if not(noref) 
	 [C,href]=contour(a,b,z',[ref ref],'LineStyle','-','LineColor',col,'LineWidth',lwid*1.5);
	end
end


% Add mask
if any(~isnan(mask(:)))
 disp('Plotting mask')
 ridgepack_mask(a,b,mask,0.5*[1 1 1]);
end

% if units of x and y are identical, then make
% the aspect ratio 1:1.
if (not(isfield(nc.(X),'units')) & not(isfield(nc.(Y),'units'))) || ...
   (isfield(nc.(X),'units') & isfield(nc.(Y),'units') && ...
    strcmp(char(nc.(X).units),char(nc.(Y).units)))
 axis equal;
 axis tight;
else
 axis normal;
end
axis xy;
box on;


% label the axes
ridgepack_labelxy(nc,X,Y);

% add title
if noref
 ridgepack_title(nc,['Contour: ',char(nc.(Z).long_name),'   (',units,')'],1)
else
 ridgepack_title(nc,['Contour: ',char(nc.(Z).long_name),'   (reference = ',num2str(ref),' ',units,')'],1)
end

hold off

% return contour handle for legends if requested
if nargout==0
 clear hpos hneg href
end

drawnow

if debug; disp(['...Leaving ',mfilename]); end

