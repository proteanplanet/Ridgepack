function ridgepack_bathymetrym

% ridgepack_bathymetrym - Adds ETOPO2 bathymetry shading to a map
%
% function ridgepack_bathymetrym
%
% Adds ETOPO2 to a map and masks out land. No inputs required.
%
% Ridgepack Version 1.0
% Andrew Roberts, Naval Postgraduate School, March 2018 (afrobert@nps.edu)


global debug;
if debug; disp(['Entering ',mfilename,'...']); end

% get current map axes
hmap=get(gcf,'CurrentAxes');
if ~ismap(hmap)
 error('Current axes must be a map')
end

% edge, grid and label color
gridcolor=0.2*[1 1 1];

% set default lats and longs
maplatlimit=[-90 90];
maplonlimit=[0 360];

% Find lat and long limits
h=gcm;
if isfield(h,'maplatlimit'); maplatlimit=sort(h.maplatlimit); end
if isfield(h,'maplonlimit'); maplonlimit=h.maplonlimit; end

% get ETOPO data file
if isempty(getenv('ETOPOFILE'))
 disp('You must set the ETOPOFILE environment variable to the local location of')
 disp('your ETOPO2.raw.bin file. You can download ETOPO2.raw.bin from the internet.')
 disp('Then set the location of ETOPOFILE in your startup.m using setenv.')
 error('You must set the ETOPOFILE environment variable to the location of ETOPO2.raw.bin')
else
 try
  [Z,refvec] = etopo(getenv('ETOPOFILE'),5,maplatlimit,maplonlimit);
 catch
  error(['Unable to find ETOPO dataset ',getenv('ETOPOFILE')]);
 end
end

% make the minimum contour less than four standard deviations 
% from the maximum in 1000's of meters
dgrad=1000;
mincont=-dgrad*(floor(max(min(Z(:)),4*std(Z(Z<0))/dgrad)));
contval=[-8000:500:-1000 -500 -250 -100 0 500];
contchoice=contval(contval>=mincont);

% generate colormap
cmap=ridgepack_colormap(contchoice,0,'bluered');

% make land grey
cmap(end,:)=0.90*[1 1 1];
colormap(cmap);

% Generate true color array
[zindex,truecolor]=ridgepack_colorindex(Z,contchoice,0);

% plot the data
hh=geoshow(gca,truecolor,refvec,'DisplayType','texturemap');

% Remove data from legend
hCGroup=hggroup;
set(hh,'Parent',hCGroup)
set(get(get(hCGroup,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');

land = shaperead('landareas.shp', 'UseGeoCoords', true);
hl=linem([land.Lat],[land.Lon],'Color',0.5*[1 1 1]);

% set colors
ridgepack_colorbar(contchoice,'meters','linear','vertical',0);

drawnow;

if debug; disp(['...Leaving ',mfilename]); end


