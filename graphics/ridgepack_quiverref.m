function ridgepack_quiverref(x,y,u,v,units,reftype,veccol,cont)

% ridgepack_quiverref - Vector plotting function for either map or Cartesian axes.
%
% function ridgepack_quiverref(x,y,u,v,units,reftype,veccol,cont)
%
% This function is a substitute for the standard version of quiver
% and quiverm available using a vanilla release of matlab. This version
% assumes a 2D vector field being plotted using a gridded flow field 
% from numerical models with a regular geometry. The function enables
% the scaling of vectors according to a reference vector plotted in the 
% lower right hand corner of the plot axes. The function works for both
% map and Cartesian axes and allows the color of vectors to be changed.
% 
% If a reference value is not provided, the reference value is calculated 
% by rounding the median, mean or maximum lengths of the quiver vectors to the 
% first significant digit. Scaling of vectors still occurs even 
% if the reference vector plotting is switched off. This enables
% different subplots to share identical scaling so that the relative
% magnitude of vectors can be compared between subplots (provided they 
% share the same grid). 
%
% The function includes the ability to plot color vectors, all of equal
% length but color coded according to their magnitude. In this case, a 
% colorbar is provided, complete with units, and no scaling vector is 
% plotted.
%
% Input:
% x       - x-coordinate or latitude
%
% y       - y-coordinate or longitude
%
% u       - u-component (Cartesian +x-direction, map +longitude-direction)
%
% v       - v-component (Cartesian +y-direction, map +latitude-direction)
%
% units   - a string providing the units of the vector field.  This assumes
%           the default Tex interpreter is being used as set 
%           elsewhere using set(0,'DefaultTextInterpreter','tex'). Optional.
%
% reftype - character variable specifying the type of reference vector.
%           Allowable values are 'median' for giving a reference
%           vector based on max(median,mean), or 'max' for giving the
%           reference vector based on the maximum.  This argument may 
%           be omitted with 'max' as the default. If reftype is entered 
%           as a number, then this is the value of the reference vector 
%           in the units of the data in u and v. If veccol is set to
%           'col', this argument has no effect and should be entered
%           only as a dummy argument. Optional.
%
% veccol  - color of the vectors to be plotted.  This may either be in the
%           form of RGB or as a single letter, such as 'b' for blue
%           using standard Matlab color specifications. This may be set
%           to 'col' if the vectors are to be color coded by magnitude
%           instead of sized by magnitude. In this case, all vectors
%           have the same size based on the optimum value for the grid
%           provided, and a colorbar is provided to reference the values.
%           Optional.
%
% cont    - contour levels to be used if color shading the vectors.  For
%           this to work, veccol='col', and as a result all vectors
%           are made of equal length but color coded according to their
%           magnitude. The contour levels must be non-zero, and must include
%           at least one value.   The values at the divider between
%           vector color categories.  Optional.
%
%
% Output:
% Output is graphical to the current active figure axes. All vectors
% are centered on the grid points, rather than the starting point of the 
% vector being positioned on the grid point. This can easily be changed
% in the code if required, but is the preference of the author.
%
% Andrew Roberts, Naval Postgraduate School, March 2018  (afrobert@nps.edu)

global debug;
if debug; disp(['Entering ',mfilename,'...']); end

% make sure the mapping toolbox is present
h=ver('map') ; if length(h)==0 ; error('Mapping toolbox not installed') ; end

% error checking of inputs
if nargin<4; error('Missing vector field input data'); end

% set default units
if nargin<5 ; units=''; end

% set default reftype value
if nargin<6 ; reftype='max' ; end

% default vector color is blue
col=false;
if nargin<7 ; 
 veccol='b' ; 
elseif ischar(veccol) & strcmp(veccol,'col')
 col=true ;
 refvec=false;
 if nargin<8; error('No contour values provided'); end
end

% get current axis 
h=get(gcf,'CurrentAxes');

% use meshgrid if needed
sx=size(x); sy=size(y);
if (sx(1)==1 & sy(1)==1) | (sx(2)==1 & sy(2)==1);
 [x,y]=meshgrid(x,y); x=x'; y=y';
elseif sx(1)==1 | sx(2)==1;
 error('Dimensions of x and y are inconsistent')
elseif sy(1)==1 | sy(2)==1;
 error('Dimensions of x and y are inconsistent')
elseif sx~=sy
 error('Dimensions of x and y are inconsistent')
end

% check sizes all agree in input data
if size(x,1)~=size(y,1) | size(x,2)~=size(y,2); 
 size(x)
 size(y)
 error('x and y sizes disagree'); 
end
if size(x,1)~=size(u,1) | size(x,2)~=size(u,2) ; 
 size(x)
 size(u)
 error('x and u sizes disagree'); 
end
if size(y,1)~=size(v,1) | size(y,2)~=size(v,2); 
 size(y)
 size(v)
 error('y and v sizes disagree'); 
end

% If plotting on a matlab map, determine if the axes are map or cartesian
% coordinates, and if the former calculate mapping to plot axis, and 
% then do vector field otherwise just plot the vector field.
if ismap(h)

 disp('Quivering on current map axes')

 % check this function can be used on the chosen map
 mapobj=gcm;
 if strmatch(mapobj.mapprojection,'globe')
   error('This function does not work with the globe projection')
 end

 % set lat and lon
 lat=x;
 lon=y;

 % get x and y location on the map
 sz=size(x);
 mstruct=gcm;
 [x,y] = mfwdtran(mstruct,lat,lon,h,'none');
 xz=size(x);
 if sz~=xz
  error('Change in size of x using mfwdtran. Try changing surface to none in the code')
 end

 % get angle on the map, but do not distort the length according to the projection
 % so that all vectors can use the same reference vector.  DO NOT project
 % the length of the vector to be different in x and y directions.
 [th,z] = cart2pol(u,v);
 [thproj,len] = vfwdtran(mstruct,lat,lon,90*ones(size(lat)));
 [u,v] = pol2cart(th+deg2rad(thproj),z);

elseif isempty(h)

  disp('Creating new Cartesian axes')

  % get the magnitude of the vector field
  [th,z] = cart2pol(u,v);

  % set up axes
  h=axes('xlim',[min(x(:)) max(x(:))],'ylim',[min(y(:)) max(y(:))]);
  axis xy;
  axis tight;
  set(gca,'Layer','bottom');
  box on;

else

 disp('Quivering on current Cartesian axes')

 % get the magnitude of the vector field
 [th,z] = cart2pol(u,v);

end

% set font size
figout=get(h,'OuterPosition');
fontsize=min(11,max(8,10*sqrt((figout(3).^2)+(figout(4).^2))/sqrt(2)));
set(h,'fontsize',fontsize);

% get current axis again
h=get(gcf,'CurrentAxes');

% turn on TeX interpreter if it is not already on
set(gcf,'DefaultTextInterpreter','tex')

% remove masked grid points and zero values from the input by filling coordinates 
% with NaN keep corner points to make sure entire grid is plotted.
xx=x;
yy=y;
x(isnan(u))=NaN;
y(isnan(u))=NaN;
x(sqrt(u.^2+v.^2)<eps)=NaN;
y(sqrt(u.^2+v.^2)<eps)=NaN;
x(1,1)=xx(1,1); y(1,1)=yy(1,1);
x(end,1)=xx(end,1); y(end,1)=yy(end,1);
x(1,end)=xx(1,end); y(1,end)=yy(1,end);
x(end,end)=xx(end,end); y(end,end)=yy(end,end);


% Scale the vectors according to the reference arrow vector length based on
% the mean distance between grid points. This is a good measure, as it remains 
% constant for multiple plots using the same grid with different values.
x1=abs(diff(x')); x2=abs(diff(x)); 
y1=abs(diff(y')); y2=abs(diff(y));
[th,z1] = cart2pol(x1,y1); [th,z2] = cart2pol(x2,y2);
scalelength=min(mean(z1(~isnan(z1))),mean(z2(~isnan(z2))));

% Calculate reference vector length based on rounded median
% or maximum value of plot.  The default is median based.
if isnumeric(reftype) & ~col
	disp('Calculating reference vector based on input number');
	refval=abs(reftype);
elseif strcmp(reftype,'median') & ~col
	disp('Calculating reference vector based on median');
        z(z==0)=NaN;
	refval=max(mean(z(~isnan(z))),median(z(~isnan(z))));
elseif strcmp(reftype,'max') & ~col
	disp('Calculating reference vector based on maximum');
	refval=max(z(~isnan(z)));
elseif ~col
	error('reftype must be either "max" or "median"');
else
	disp('Color vectors being used with constant length');
end

% Remove NaN values that will not be plotted
% and turn points into a row of coordinates
u=u(~isnan(x))';
v=v(~isnan(x))';
y=y(~isnan(x))';
x=x(~isnan(x))';

% Set arrow size (1 = full length of vector)
arrow=0.40;

% set the length of the vectors to be constant and generate 
% a color for each particular vector
if col 
 
 % set line width of vector
 vecwidth=0.4;

 % Set accenture of arrows (1 = full arrow)
 %accenture=0.30;
 accenture=1.0;

 % get the magnitude and direction of each vector
 [th,z] = cart2pol(u,v);

 % check cont has at least one value and is non-zero
 cont=cont(cont>0);
 if isempty(cont) | length(cont)==0 | any(cont==0)
  error('cont must be non-zero and contain at least one value')
 end

 % arrange contour values from min to max, and 
 % add an extra value for the purpose of processing
 cont=sort(abs(cont));
 cont(end+1)=cont(end)+1;

 % set the colormap 
 hc=colormap(jet(length(cont)));

 for i=1:length(cont)

  if i==1
   mask=find(z<cont(i));
  elseif i==length(cont)
   mask=find(z>=cont(i-1));
  else
   mask=find(z<cont(i) & z>=cont(i-1));
  end

  % Center vectors over grid points
  [u,v] = pol2cart(th(mask),scalelength);
  xstart=x(mask)-0.5*u;
  xend=x(mask)+0.5*u;
  ystart=y(mask)-0.5*v;
  yend=y(mask)+0.5*v;

  % Get x coordinates of each vector plotted
  lx = [xstart; ...
       xstart+(1-arrow*accenture)*(xend-xstart); ...
       xend-arrow*(u+arrow*v); ...
       xend; ...
       xend-arrow*(u-arrow*v); ...
       xstart+(1-arrow*accenture)*(xend-xstart); ...
       xstart];

  % Get y coordinates of each vector plotted
  ly = [ystart; ...
       ystart+(1-arrow*accenture)*(yend-ystart); ...
       yend-arrow*(v-arrow*u); ...
       yend; ...
       yend-arrow*(v+arrow*u); ...
       ystart+(1-arrow*accenture)*(yend-ystart); ...
       ystart]; 

  % set height just above land outline in map
  if ismap(h)
   lz=0.00002*ones(size(ly));
   lz(isnan(ly))=NaN;
  else
   lz=zeros(size(ly));
   lz(isnan(ly))=NaN;
  end

  % patch the vectors
  patch(lx,ly,lz,squeeze(hc(i,:)),'EdgeColor',squeeze(hc(i,:)),'LineWidth',vecwidth);

 end

 % generate the colorbar 
 hcb=ridgepack_colorbar([0 cont],units,'linear','vertical',0);

else 

 % set line width of vector
 vecwidth=0.4;

 % Set accenture of arrows (1 = full arrow)
 accenture=0.0;

 % set scale value based on refval and scale length
 roundp=floor(log10(refval));
 refval=floor(refval/(10^roundp))*(10^roundp);
 scale=scalelength/refval;

 % Ensure that the data aspect ratio is 1:1 in the x:y direction
 set(h,'DataAspectRatio',[1 1 1]);

 % Center vectors over grid points
 xstart=x-0.5*scale*u;
 xend=x+0.5*scale*u;
 ystart=y-0.5*scale*v;%*ratio;
 yend=y+0.5*scale*v;%*ratio;


 % Get x coordinates of each vector plotted
 lx = [xstart; x; ...
      xstart+(1-arrow*accenture)*(xend-xstart); ...
      xend-arrow*(scale*u+arrow*(scale*v)); ...
      xend; ...
      xend-arrow*(scale*u-arrow*(scale*v)); ...
      xstart+(1-arrow*accenture)*(xend-xstart); ...
      repmat(NaN,size(x))];

 % Get y coordinates of each vector plotted
 ly = [ystart; y; ...
      ystart+(1-arrow*accenture)*(yend-ystart); ...
      yend-arrow*(scale*v-arrow*(scale*u)); ...
      yend; ...
      yend-arrow*(scale*v+arrow*(scale*u)); ...
      ystart+(1-arrow*accenture)*(yend-ystart); ...
      repmat(NaN,size(y))];

 % set height just above land outline in map
 if ismap(h)
  lz=0.00002*ones(size(ly));
  lz(isnan(ly))=NaN;
 else
  lz=zeros(size(ly));
  lz(isnan(ly))=NaN;
 end

 % Plot the vectors
 if debug
  plot(x,y,'r.')
 else
  line(lx,ly,lz,'Color',veccol,'LineWidth',vecwidth);
 end

 % Draw the reference vector key at altitude 2 above the map and grid
 ridgepack_vecref(gca,scalelength,arrow,refval,units,veccol,vecwidth,accenture,lx,ly);

end

if debug; disp(['...Leaving ',mfilename]); end

