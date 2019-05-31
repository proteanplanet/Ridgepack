function [h]=ridgepack_stereom(hemisphere,grid,label)

% ridgepack_stereom - Plots a native polar stereographic map
%
% function [h]=ridgepack_stereom(hemisphere,grid,label)
%
% This function plots a native polar stereographic map for the Arctic
% and Antarctic that does not use the MATLAB Mapping Toolbox. As a 
% result, it may be used with a basic MATLAB version. 
%
% INPUT:
%
% hemisphere - set to 1 for the Arctic, and -1 for the Antarctic (default=1)
% grid       - logical switching on or off plotting of a grid (default=true)
% label      - logical switching on or off labeling of the grid (default=true)
%
% OUTPUT:
%
% h - handle for the current axes
%
% Ridgepack Version 1.1
% Andrew Roberts, Los Alamos National Laboratory, May 2019 (afroberts@lanl.gov)

global debug;
if debug; disp(['Entering ',mfilename,'...']); end

% error checking
if nargin>=1 & hemisphere~=1 & hemisphere~=-1
 error('hemisphere must be either 1 or -1')
elseif nargin==0
 hemisphere=1;
end

if nargin>=2 & ~islogical(grid)
 error('grid must be true or false')
elseif nargin<=1
 grid=true;
end

if nargin>=3 & ~islogical(label)
 error('label must be true or false')
elseif nargin<=2
 label=true;
end

% generate SSM/I grid
if hemisphere==1
 nc=ridgepack_gridgen('',8);
elseif hemisphere==-1
 nc=ridgepack_gridgen('',9);
else
 error('hemisphere specificiation error')
end
[x,y]=meshgrid(nc.x.data,nc.y.data);

plot(x(:,1),y(:,1),'k');
hold on
plot(x(:,end),y(:,end),'k');
plot(x(1,:),y(1,:),'k');
hold on
plot(x(:,1),y(:,1),'k');
plot(x(:,end),y(:,end),'k');

xmins=min(x(:));
xmaxs=max(x(:));
ymins=min(y(:));
ymaxs=max(y(:));

clear nc x y

axis tight
axis equal
set(gca,'XTick',[],'YTick',[])

if grid

 for llat=hemisphere*[0:10:80]

  lon=[0:0.01:360];
  lat=llat*ones(size(lon));

  [X,Y]=ridgepack_geodetictoxy(lat,lon,1);

  plot3(X,Y,0.01*ones(size(X)),'Color',0.75*[1 1 1],'LineStyle',':')

  if hemisphere==1 & llat>40 & llat<80 & label
   lon=[135];
   lat=llat*ones(size(lon));
   for ll=1:length(lon)
    [X,Y]=ridgepack_geodetictoxy(lat(ll),lon(ll),1);
    text(X,Y,0.01,[num2str(lat(ll)),'$^{\circ}$N'],...
        'FontSize',7,'Color',0.75*[1 1 1],...
        'HorizontalAlignment','center',...
        'VerticalAlignment','bottom',...
        'Rotation',lon(ll)+180);
   end
  elseif hemisphere==-1 & llat<-40 & llat>-80 & label
   lon=[50];
   lat=llat*ones(size(lon));
   for ll=1:length(lon)
    [X,Y]=ridgepack_geodetictoxy(lat(ll),lon(ll),1);
    text(X,Y,0.01,[num2str(abs(lat(ll))),'$^{\circ}$N'],...
        'FontSize',7,'Color',0.75*[1 1 1],...
        'HorizontalAlignment','center',...
        'VerticalAlignment','bottom',...
        'Rotation',lon(ll)-90-10);
   end
  end
 
 end

 for llon=0:30:330

  if llon==0 | llon==90 | llon==180 | llon==270
   lat=hemisphere*[0:0.01:90];
  else
   lat=hemisphere*[0:0.01:80];
  end
  lon=llon*ones(size(lat));

  [X,Y]=ridgepack_geodetictoxy(lat,lon,1);

  plot3(X,Y,0.01*ones(size(X)),'Color',0.75*[1 1 1],'LineStyle',':')

  if hemisphere==1 & (llon==0 | llon==180) & label 
   lon=llon;
   if llon==0
    lat=[46];
   elseif llon==180
    lat=[43];
   end
   for ll=1:length(lon)
    [X,Y]=ridgepack_geodetictoxy(lat(ll),lon(ll),1);
    text(X,Y,0.01,[num2str(lon(ll)),'$^{\circ}$E'],...
        'FontSize',7,'Color',0.75*[1 1 1],...
        'HorizontalAlignment','center',...
        'VerticalAlignment','bottom',...
        'Rotation',-90);
   end
  elseif hemisphere==-1 & (llon==0 | llon==180) & label 
   lon=llon;
   if llon==0
    lat=[-57];
   elseif llon==180
    lat=[-55];
   end
   for ll=1:length(lon)
    [X,Y]=ridgepack_geodetictoxy(lat(ll),lon(ll),1);
    text(X,Y,0.01,[num2str(lon(ll)),'$^{\circ}$E'],...
        'FontSize',7,'Color',0.75*[1 1 1],...
        'HorizontalAlignment','center',...
        'VerticalAlignment','bottom',...
        'Rotation',-90);
   end
  end

 end

end

%landcheck=true;
landcheck=false;
if landcheck

 gridloclr='/Volumes/CICESatArray/E3SM/lrhrequiv/grid';
 cd(gridloclr)

 landm = shaperead('E3SM_V1_C_grid_Coast.shp','UseGeoCoords',true);
 if hemisphere==1
  lat=landm.Lat(landm.Lat(:)>0 | isnan(landm.Lat(:)));
  lon=landm.Lon(landm.Lat(:)>0 | isnan(landm.Lat(:)));
 elseif hemisphere==-1
  lat=landm.Lat(landm.Lat(:)<0 | isnan(landm.Lat(:)));
  lon=landm.Lon(landm.Lat(:)<0 | isnan(landm.Lat(:)));
 end
 [X,Y]=ridgepack_geodetictoxy(lat,lon,1);
 h1=plot3(X,Y,0.005*ones(size(X)),'Color',[0.5 0.5 0.5]);

end

if hemisphere==-1
 axis ij
end

% set plot axis limits
xlim([xmins xmaxs]);
ylim([ymins ymaxs]);

% Get handle
h=gca;


