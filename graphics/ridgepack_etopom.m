function [nc]=ridgepack_etopom(sample,maplatlimit)

% ridgepack_etopom - Adds ETOPO2 relief shading to a map in color or output to netcdf
% 
% function [nc]=ridgepack_etopom(sample,maplatlimit)
%
% This adds ETOPO2 data to a map, in color, sampled at every 'sample' data 
% point. A map must first be drawn before running this function. If nc is 
% specified as output, an nc structure containing the sampled etopo data is 
% created instead of plotting the data, in which case a map need not be drawn.
%
% INPUT:
%
% Sample - sample rate of 2 minute data.  The default is every fifth point.
% maplatlimit - a two element vector with the min and max latitudes to be extracted.
%          This is only read if an nc structure is being created, otherwise
%          the information is automatically obtained from the map to which
%          etopo data is being added.
%
%
% OUTPUT:
%
% nc - etopo data on (lat x long) grid.  If output is requested, then
%     the data is not plotted graphically.
%
% Ridgepack Version 1.0
% Andrew Roberts, Naval Postgraduate School, March 2018 (afrobert@nps.edu)


global debug;
if debug; disp(['Entering ',mfilename,'...']); end

% check inputs and set defaults
if nargin==0; 
 sample=5; 
elseif ~isnumeric(sample) || length(sample)~=1
 sample
 error('Sample must be a number')
end

% Find lat and long limits
if nargin<=1;
 maplatlimit=[-90 90];
elseif ~isnumeric(maplatlimit) || length(maplatlimit)~=2
 maplatlimit
 error('maplatlimit is incorrect')
end
maplonlimit=[0 360];

% Get map handle
if nargout==0 && ~ismap(gca)
 error('A map must first be drawn')
elseif nargout==0 && nargin<=1;
 h=gcm;
 if isfield(h,'maplatlimit'); maplatlimit=sort(h.maplatlimit); end
 if isfield(h,'maplonlimit'); maplonlimit=h.maplonlimit; end
end

% Edge, grid and label color
gridcolor=0.2*[1 1 1] ;

% Set default color gradation
dgrad=500;

% get ETOPO data file
if isempty(getenv('ETOPOFILE'))
 disp('You must set the ETOPOFILE environment variable to the local location of')
 disp('your ETOPO2.raw.bin file. You can download ETOPO2.raw.bin from the internet.')
 disp('Then set the location of ETOPOFILE in your startup.m using setenv.')
 error('You must set the ETOPOFILE environment variable to the location of ETOPO2.raw.bin')
else
 try
  [Z,refvec] = etopo(getenv('ETOPOFILE'),sample,maplatlimit,maplonlimit);
 catch
  error(['Unable to find ETOPO dataset ',getenv('ETOPOFILE')]);
 end
end

if nargout==0;  % plot data on a map

 % get current map axes
 hmap=get(gcf,'CurrentAxes');
 if ~ismap(hmap)
  error('Current axes must be a map')
 end

 % find min and max contour 
 mincont=-dgrad*(floor(max(min(Z(:)),4*std(Z(Z<0))/dgrad))); 
 maxcont=dgrad*(floor(min(max(Z(:)),4*std(Z(Z>0))/dgrad))); 
 contval=[-8000:1000:-1000 -500 -250 0 250 500 1000:1000:8000];
 %contchoice=contval(contval>=mincont & contval<=maxcont);
 contchoice=[-6000:1000:-1000 -500 -250 0 250 500 1000:1000:8000];

 % generate colormap
 cmap=ridgepack_colormap(contchoice,0,'bluered');

 % Generate true color array
 [zindex,truecolor]=ridgepack_colorindex(Z,contchoice,0);

 % plot the data
 hh=geoshow(gca,truecolor,refvec,'DisplayType','texturemap');

 % Remove from legend
 hCGroup=hggroup;
 set(hh,'Parent',hCGroup)
 set(get(get(hCGroup,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');

 % add in coastline
 land = shaperead('landareas.shp', 'UseGeoCoords', true);
 hl=linem([land.Lat],[land.Lon],'Color',0.5*[1 1 1]);

 % set colors
 ridgepack_colorbar(contchoice,'meters','linear','vertical',0);

else % write netcdf file

 xx=sample*2;
 nc.attributes.title='etopo2 subsampled on 10 minute grid';
 nc.latitude.data=min(maplatlimit(:))+(xx/120):(xx/60):max(maplatlimit(:))-(xx/120);
 nc.latitude.long_name='latitude';
 nc.latitude.units='degrees_north';
 nc.latitude.dimension={'latitude'};
 nc.longitude.data=-180+(xx/120):(xx/60):180-(xx/120);
 nc.longitude.long_name='longitude';
 nc.longitude.units='degrees_east';
 nc.longitude.dimension={'longitude'}';
 nc.Z.data=Z;
 nc.Z.long_name='etopo';
 nc.Z.units='m';
 nc.Z.dimension={'latitude','longitude'};
 nc=ridgepack_struct(nc); 

end

drawnow

if debug; disp(['...Leaving ',mfilename]); end


