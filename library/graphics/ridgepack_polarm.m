function [maph]=ridgepack_polarm(varargin)

% ridgepack_polarm - Square or rectangular polar stereographic map with land and grid
%
% function [maph]=ridgepack_polarm(vargin)
%
% This function generates a polar stereographic for either the 
% Arctic [default] or Antarctic on a square-framed map (the default
% stereographic map in matlab's mapping toolbox produces a circular plot,
% but this is not always the most space-efficient method for 
% displaying oceanographic and meteorological information at the 
% poles). The are also several special Arctic preset maps available in 
% this function for displaying information about the Arctic System.
% 
% Several optional vargin arguments may be used:
%
% INPUT:
%
% 'lat'   - Provides the maximum equatorward extent the map
%           and is proceeded by a scalar latitude.  Positive
%           latitudes produce an Arctic polar stereographic
%           map and negative latitudes produce an Antarctic 
%           polar stereographic map.  The default is 60 degrees
%           north, giving a default Arctic map.
%
%  lat    - Actual latitude value following 'lat'.
%
% 'lon'   - Central meridian in the lower (upper) portion for the
%           northern (southern) hemisphere. The default is the
%           Greenwich Meridian.
%
%  lon    - Actual value of lon following 'lon'.
%
% 'noland'- Stop land mask being plotted.
%
% 'grid'  - Overlays a grid on the map.
%
% 'label' - Adds latitude and longitude labels to the grid.
%
% 'asm'   - Arctic System map projection designed to include
%           areas that may be considered part of the regional
%           "Arctic System". This option does not work with 'seaice'. 
%
% 'rasm'  - Map of the entire domain of the Regional Arctic System Model
%
% 'rasmshort' - Nearly entire domain of the Regional Arctic System Model
%
% 'rasmlong' - Encompass entire domain of the Regional Arctic System Model
%
% 'seaice'- Plot most efficient rectangular map for Arctic
%           sea ice. This option does not work if 'asm' is selected.
%
% 'seaicerasm'- Plot most efficient rectangular map for Arctic
%           sea ice in rasm orientation. Does not work if 'asm' is selected.
%
% 'centralarctic' - Plot of central arctic with Queen Elizabeth Islands
%           positioned parallel to the bottom of the square map frame. This
%           map is most relevant for analyzing arctic sea ice mechanics.
%
% 'centralarctic2' - A zoomed out version of centralarctic.
%
% 'rasmmask' - Plot of central arctic that includes the central arctic RASM mask.
%
% 'shipping' - Map of central Arctic applicable to shipping
%
% 'antarctic' - Map for Antarctic sea ice zone.
%
%
% OUTPUT:
%
% This function generates a handle for the map:  maph
%
% An Azimuthal Stereographic Projection map projection is
% used.  Typing 'help stereo' in matlab will provide more
% information about this projection.  To see the distortion 
% provided by this map, use the tissot function.
%
% Examples:
% ridgepack_polarm('asm') 
% ridgepack_polarm('seaice')
% ridgepack_polarm('lat',40,'lon',-40,'grid','label')
% ridgepack_polarm('lat',-60,'grid')
%
% When printing maps with Postscript or EPSC, it is advisable to use the 
% -zbuffer option to avoid thin white filaments showing up as part of the land 
% filling. This is an artifact from within Matlab's Mapping toolbox.
%
% Ridgepack Version 1.0
% Andrew Roberts, Naval Postgraduate School, March 2018 (afrobert@nps.edu)


global debug;
if debug; disp(['Entering ',mfilename,'...']); end

% Check that the mapping toolbox is installed
h=ver('map') ; 
if length(h)==0 ; error('Mapping toolbox not installed') ; end

% defaults
centralmeridian=0 ; 
equatorextent=60 ; 
parallelstop=80;
grid=0;
label=0;
land=1;
asm=0;
rasm=0;
rasmshort=0;
rasmlong=0;
shipping=0;
centralarctic=0;
centralarctic2=0;
rasmmask=0;
seaice=0;
seaicerasm=0;
antarctic=0;

% edge, grid and label color
gridcolor=0.6*[1 1 1] ; 

if nargin == 0
 disp('Default Arctic stereographic map');
else
 i=0;
 while i < nargin
  i=i+1;
  switch lower(varargin{i})
   case 'seaice'
	   seaice=1;
   case 'seaicerasm'
	   seaicerasm=1;
   case 'shipping'
	   shipping=1;
   case 'asm'
	   asm=1;
   case 'rasm'
	   rasm=1;
   case 'rasmshort'
	   rasmshort=1;
   case 'rasmlong'
	   rasmlong=1;
   case 'centralarctic'
	   centralarctic=1;
   case 'centralarctic2'
	   centralarctic2=1;
   case 'rasmmask'
	   rasmmask=1;
   case 'antarctic'
	   antarctic=1;
   case 'lat'
  	   if(i+1>nargin || ischar(varargin{i+1}))
	    error(['Missing latitude argument for ',varargin{i}])
  	   end
	   i=i+1;
	   equatorextent=varargin{i};
   case 'lon'
  	   if(i+1>nargin || ischar(varargin{i+1}))
	    error(['Missing longitiude argument for ',varargin{i}])
  	   end
	   i=i+1;
	   centralmeridian=varargin{i};
   case 'grid'
	   grid=1;
   case 'label'
	   label=1;
   case 'noland'
	   land=0;
   otherwise
	   error(['Option ',varargin{i},' is incorrect'])
  end
 end
end

if asm % create preset arctic system map
 equatorextent=40;
 centralmeridian=-45;
elseif seaice % create preset arctic system map
 equatorextent=40;
 centralmeridian=-45;
elseif seaicerasm
 equatorextent=40;
 centralmeridian=-114;
elseif shipping
 equatorextent=51;
 centralmeridian=-114;
elseif rasm
 equatorextent=30;
 centralmeridian=-114;
elseif rasmshort
 equatorextent=35;
 centralmeridian=-114;
elseif rasmlong
 equatorextent=10;
 centralmeridian=-114;
elseif centralarctic
 equatorextent=60;
 centralmeridian=-58;
elseif centralarctic2
 equatorextent=60;
 centralmeridian=-58;
elseif rasmmask
 equatorextent=60;
 centralmeridian=-58;
elseif antarctic
 equatorextent=-55;
 centralmeridian=0;
end

% make polar map southern hemisphere if southern extent is zero
if(abs(equatorextent)>89)
	error('Latitude extent of map too close to a pole.');
elseif(equatorextent<0)
	mult=-1;
	if debug; disp('Southern Hemisphere'); end
else
	mult=1;
	if debug; disp('Northern Hemisphere'); end
end

% clear current axes
%cla; 

% set default map properties
defaultm; 

% map projection
maph=axesm('MapProjection','stereo',...
  'AngleUnits','degrees',...
  'Aspect','normal',...
  'FalseNorthing',0,...
  'FalseEasting',0,...
  'MapLatLimit',sort([0 mult*90]),...
  'Geoid',[1 0],...
  'Origin',[mult*90 centralmeridian 0],...
  'Scalefactor',1,...
  'Frame','on',...
  'FFill',2000,...
  'FLatLimit',[-mult*Inf mult*90],...
  'FEdgeColor','none',...
  'FFaceColor','none',...
  'FLineWidth',0.25);

% fit map tightly to axes
tightmap 
xlim=get(gca,'xlim');
ylim=get(gca,'ylim');
[lat,leftlon]=minvtran(xlim(1),ylim(1)+0.5*diff(ylim));
[lat,rightlon]=minvtran(xlim(2),ylim(1)+0.5*diff(ylim));
[lat,bottomlon]=minvtran(xlim(1)+0.5*diff(xlim),ylim(1));
[lat,toplon]=minvtran(xlim(1)+0.5*diff(xlim),ylim(2));
if asm
 [x1,y] = mfwdtran(equatorextent+7,leftlon);
 [x2,y] = mfwdtran(equatorextent+9,rightlon);
 [x,y1] = mfwdtran(equatorextent+3,bottomlon);
 [x,y2] = mfwdtran(equatorextent,toplon);
elseif seaice
 [x1,y] = mfwdtran(equatorextent+16,leftlon);
 [x2,y] = mfwdtran(equatorextent+17,rightlon);
 [x,y1] = mfwdtran(equatorextent+5,bottomlon);
 [x,y2] = mfwdtran(equatorextent-0,toplon);
elseif seaicerasm
 [x1,y] = mfwdtran(equatorextent-1,leftlon);
 [x2,y] = mfwdtran(equatorextent+10,rightlon);
 [x,y1] = mfwdtran(equatorextent+16,bottomlon);
 [x,y2] = mfwdtran(equatorextent+23,toplon);
elseif shipping
 [x1,y] = mfwdtran(equatorextent+12,leftlon);
 [x2,y] = mfwdtran(equatorextent+10,rightlon);
 [x,y1] = mfwdtran(equatorextent+16,bottomlon);
 [x,y2] = mfwdtran(equatorextent+17,toplon);
elseif rasm
 [x1,y] = mfwdtran(equatorextent-4.5,leftlon);
 [x2,y] = mfwdtran(equatorextent+10.3,rightlon);
 [x,y1] = mfwdtran(equatorextent+15,bottomlon);
 [x,y2] = mfwdtran(equatorextent+16.1,toplon);
elseif rasmshort
 [x1,y] = mfwdtran(equatorextent-5,leftlon);
 [x2,y] = mfwdtran(equatorextent+10,rightlon);
 [x,y1] = mfwdtran(equatorextent+13,bottomlon);
 [x,y2] = mfwdtran(equatorextent+13,toplon);
elseif rasmlong
 [x1,y] = mfwdtran(equatorextent+7.5,leftlon);
 [x2,y] = mfwdtran(equatorextent+22,rightlon);
 [x,y1] = mfwdtran(equatorextent+19,bottomlon);
 [x,y2] = mfwdtran(equatorextent+20,toplon);
elseif centralarctic
 [x1,y] = mfwdtran(equatorextent+9,leftlon);
 [x2,y] = mfwdtran(equatorextent+23,rightlon);
 [x,y1] = mfwdtran(equatorextent+22.5,bottomlon);
 [x,y2] = mfwdtran(equatorextent+13,toplon);
elseif centralarctic2
 [x1,y] = mfwdtran(equatorextent+6,leftlon);
 [x2,y] = mfwdtran(equatorextent+19,rightlon);
 [x,y1] = mfwdtran(equatorextent+21,bottomlon);
 [x,y2] = mfwdtran(equatorextent+10,toplon);
elseif rasmmask
 [x1,y] = mfwdtran(equatorextent+6.5,leftlon);
 [x2,y] = mfwdtran(equatorextent+19.8,rightlon);
 [x,y1] = mfwdtran(equatorextent+22.5,bottomlon);
 [x,y2] = mfwdtran(equatorextent+11,toplon);
elseif antarctic
 [x1,y] = mfwdtran(equatorextent-7,leftlon);
 [x2,y] = mfwdtran(equatorextent-1,rightlon);
 [x,y1] = mfwdtran(equatorextent-2.5,bottomlon);
 [x,y2] = mfwdtran(equatorextent+2.5,toplon);
else
 [x1,y] = mfwdtran(equatorextent,leftlon);
 [x2,y] = mfwdtran(equatorextent,rightlon);
 [x,y1] = mfwdtran(equatorextent,bottomlon);
 [x,y2] = mfwdtran(equatorextent,toplon);
end

set(gca,'Xlim',[x1 x2]); 
set(gca,'Ylim',[y1 y2]);


[latc1,lonc1]=minvtran(x1,y1);
[latc2,lonc2]=minvtran(x1,y2);
[latc3,lonc3]=minvtran(x2,y2);
[latc4,lonc4]=minvtran(x2,y1);

% Set fontsize
figout=get(gca,'Position');
%fontsize=max(4,min(8,10*sqrt((figout(3).^2)+(figout(4).^2))));
fontsize=max(5,min(8,10*sqrt((figout(3).^2)+(figout(4).^2))));
set(gca,'fontsize',fontsize);

% Scaled degree spacing of grid dots 
ddeg=distance(latc1,lonc1,latc3,lonc3)*sqrt(2)/(360*sqrt((figout(3).^2)+(figout(4).^2)));

markersize=1;
%markersize=min(1,sqrt((figout(3).^2)+(figout(4).^2)));


set(gca,'box','on')

% add land to plot
if (land==1)
 landm = shaperead('landareas.shp', 'UseGeoCoords', true);

 % land mask 
 p1=patchm([landm.Lat],[landm.Lon],-0.00001,0.96*[1 1 1],...
           'Clipping','on','EdgeColor','none');

 % land outline (positioned at a slight altitude to remain visible)
 l1=linem([landm.Lat],[landm.Lon],0.00001,'Color',0.45*[1 1 1],'LineWidth',0.3);
end

% grid label
if(label==1)

 % meridians of latitude and longitude labels
 if centralarctic | centralarctic2 | rasmmask
  MLabelLocation=[90];
  PLabelMeridian=100;
 elseif(mult>0)
  MLabelLocation=[0 180];
  PLabelMeridian=105;
 elseif antarctic
  MLabelLocation=[0];
  PLabelMeridian=45;
 else
  MLabelLocation=[-90 90];
  PLabelMeridian=135;
 end

 % parallels to be labelled 
 if abs(equatorextent)<20
  PLabelLocation=mult*[0 30 60];
 elseif centralarctic | centralarctic2 | rasmmask
  PLabelLocation=[70];
 else
  PLabelLocation=mult*[30 40 50 60 70];
 end

 for i=1:length(MLabelLocation)
  mlon(i)=MLabelLocation(i);
  lat=mult*[5:10:85];
  clear x y xdiff1 xdiff2 ydiff1 ydiff2
  for j=1:length(lat)
   [x(j),y(j)]=mfwdtran(lat(j),mlon(i));
   xdiff1(j)=abs(x(j)-x1);
   xdiff2(j)=abs(x2-x(j));
   ydiff1(j)=abs(y(j)-y1);
   ydiff2(j)=abs(y2-y(j));
  end
  xdiff1(x<=x1 | x>=x2 | y<=y1 | y>=y2)=realmax;
  xdiff2(x<=x1 | x>=x2 | y<=y1 | y>=y2)=realmax;
  ydiff1(x<=x1 | x>=x2 | y<=y1 | y>=y2)=realmax;
  ydiff2(x<=x1 | x>=x2 | y<=y1 | y>=y2)=realmax;
  if min(xdiff1)==min([xdiff1 xdiff2 ydiff1 ydiff2])
   li=find(xdiff1==min(xdiff1));
  elseif min(xdiff2)==min([xdiff2 ydiff1 ydiff2])
   li=find(xdiff2==min(xdiff2));
  elseif min(ydiff1)==min([ydiff1 ydiff2])
   li=find(ydiff1==min(ydiff1));
  else
   li=find(ydiff2==min(ydiff2));
  end

  mlat(i)=lat(li);

  % set rotation of labels
  %if antarctic
  % rot=0;
  %else
  % rot=mod(mlon(i),180)-centralmeridian-90;
  %end
  rot=0;

  if seaicerasm
   ht=text(x(li),y(li),['$\;',num2str(mlon(i)),'^{\circ}$ '],...
      'HorizontalAlignment','center',...
      'VerticalAlignment','middle',...
      'Interpreter','latex',...
      'FontSize',fontsize,...
      'Margin',1,...
      'Rotation',rot,...
      'Color',0*[1 1 1]);
  else
   ht=text(x(li),y(li),['$\;',num2str(mlon(i)),'^{\circ}$ '],...
      'HorizontalAlignment','center',...
      'VerticalAlignment','middle',...
      'Interpreter','latex',...
      'FontSize',fontsize,...
      'Margin',1,...
      'Rotation',rot,...
      'Color',gridcolor);
  end
  for k=1:5
   extent=get(ht,'extent');
   if extent(1)<x1 | extent(2)<y1
    mlat(i)=mlat(i)+mult*10;
    [x(li),y(li)]=mfwdtran(mlat(i),mlon(i));
   elseif extent(1)+extent(3)>x2 | extent(2)+extent(4)>y2
    mlat(i)=mlat(i)+mult*10;
    [x(li),y(li)]=mfwdtran(mlat(i),mlon(i));
   end
   set(ht,'position',[x(li) y(li) 0]);
  end
  
  xminl(i)=extent(1);
  xmaxl(i)=extent(1)+extent(3);
  yminl(i)=extent(2);
  ymaxl(i)=extent(2)+extent(4);

 end
 
 Ppos=PLabelMeridian;

 plat=PLabelLocation;
 plon=Ppos*ones(size(plat));

 [x,y]=mfwdtran(plat,plon); 

 x(x<x1)=NaN;
 x(x>x2)=NaN;
 y(y<y1)=NaN;
 y(y>y2)=NaN;

 for i=1:length(x)
  if plat(i)>0
   tstring=['$',num2str(plat(i)),'^{\circ}$N'];
  elseif plat(i)<0
   tstring=['$',num2str(abs(plat(i))),'^{\circ}$S'];
  else
   tstring=['$',num2str(plat(i)),'^{\circ}$'];
  end
  if seaicerasm
   ht=text(x(i),y(i),tstring,...
      'HorizontalAlignment','center',...
      'VerticalAlignment','middle',...
      'Interpreter','latex',...
      'Margin',1,...
      'FontSize',fontsize,...
      'Color',0*[1 1 1]);
  else
   ht=text(x(i),y(i),tstring,...
      'HorizontalAlignment','center',...
      'VerticalAlignment','middle',...
      'Interpreter','latex',...
      'Margin',1,...
      'FontSize',fontsize,...
      'Color',gridcolor);
  end
  for k=1:5
   extent=get(ht,'extent');
   if extent(1)<x1 | extent(2)<y1
    mlat(i)=mlat(i)+mult*10;
    [x(i),y(i)]=mfwdtran(mlat(i),mlon(i));
   elseif extent(1)+extent(3)>x2 | extent(2)+extent(4)>y2
    mlat(i)=mlat(i)+mult*10;
    [x(i),y(i)]=mfwdtran(mlat(i),mlon(i));
   end
   set(ht,'position',[x(i) y(i) 0]);
  end
  xmink(i)=extent(1);
  xmaxk(i)=extent(1)+extent(3);
  ymink(i)=extent(2);
  ymaxk(i)=extent(2)+extent(4);
 end

else

 mlon=NaN;
 mlat=NaN;
 plon=NaN;
 plat=NaN;

end

% map grid
if(grid==1)
 for i=[0 90 -90 180]
  lat=[-90:ddeg:90];
  lon=i*ones(size(lat));
  [x,y]=mfwdtran(lat,lon);
  ii=find(i==mlon);
  if ~isempty(ii)
   idx=find(x>xminl(ii(1)) & x<xmaxl(ii(1)) & y>yminl(ii(1)) & y<ymaxl(ii(1)));
   x(idx)=NaN;
   y(idx)=NaN;
  end
  z=0.001*ones(size(x));
  if seaicerasm
   plot3(x,y,z,'k.','MarkerSize',markersize)
  else
   plot3(x,y,z,'.','Color',gridcolor,'MarkerSize',markersize)
  end
 end

 for i=[30 60 120 150 -30 -60 -120 -150]
  lat=[-80:ddeg:80];
  lon=i*ones(size(lat));
  [x,y]=mfwdtran(lat,lon);
  z=0.001*ones(size(x));
  if seaicerasm
   plot3(x,y,z,'k.','MarkerSize',markersize)
  else
   plot3(x,y,z,'.','Color',gridcolor,'MarkerSize',markersize)
  end
 end

 for j=[-80:10:80]
  lon=[-180:ddeg/cosd(j):180];
  lat=j*ones(size(lon));
  [x,y]=mfwdtran(lat,lon);
  jj=find(j==plat);
  if ~isempty(jj)
   idx=find(x>xmink(jj(1)) & x<xmaxk(jj(1)) & y>ymink(jj(1)) & y<ymaxk(jj(1)));
   x(idx)=NaN;
   y(idx)=NaN;
  end
  z=0.001*ones(size(x));
  if seaicerasm
   plot3(x,y,z,'k.','MarkerSize',markersize)
  else
   plot3(x,y,z,'.','Color',gridcolor,'MarkerSize',markersize)
  end
 end
end

%[x,y]=mfwdtran(mult*90,0);
%plot(x,y,'.','Color',gridcolor)

% Remove meridians and parallels from the legend (if one is used)
h1=handlem('parallel');
set(h1,'Clipping','on');
hCGroup=hggroup;
set(h1,'Parent',hCGroup);
set(get(get(hCGroup,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');

h2=handlem('meridian');
set(h2,'Clipping','on');
hSGroup=hggroup;
set(h2,'Parent',hSGroup);
set(get(get(hSGroup,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');

h2=handlem('patch');
set(h2,'Clipping','on');
hSGroup=hggroup;
set(h2,'Parent',hSGroup);
set(get(get(hSGroup,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');

h2=handlem('line');
set(h2,'Clipping','on');
hSGroup=hggroup;
set(h2,'Parent',hSGroup);
set(get(get(hSGroup,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');

h2=handlem('PLabel');
set(h2,'Clipping','on');
hSGroup=hggroup;
set(h2,'Parent',hSGroup);
set(get(get(hSGroup,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');

h2=handlem('MLabel');
set(h2,'Clipping','on');
hSGroup=hggroup;
set(h2,'Parent',hSGroup);
set(get(get(hSGroup,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');

h2=handlem('Frame');
set(h2,'Clipping','on');
hSGroup=hggroup;
set(h2,'Parent',hSGroup);
set(get(get(hSGroup,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');

% This line necessary to get a snug fit with horizontal colorbars
set(gca,'XTickLabel',[])
set(gca,'YTickLabel',[])


% remove map handle if not requested
if nargout==0 ; clear maph; end

z=0.002;
plotbox=plot3([x1 x2 x2 x1 x1],[y1 y1 y2 y2 y1],[z z z z z],'Color','k');
uistack(plotbox,'top')

drawnow

if debug; disp(['...Leaving ',mfilename]); end


