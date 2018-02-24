function [maph]=ridgepack_globem(varargin)

% ridgepack_globem - Generates a globe map view of earth from space
%
% function [maph]=ridgepack_globem(vlat,vlon,noland,transparent)
%
% This generates a spherical map of the earth as seen 
% from space directly above the coordinates vlat,vlon.
%
% INPUT:
%
% 'lon'       - Longitude of center of projection.
%
%  lon        - Actual lon value following 'lon'.
%
% 'lat'       - Latitude of center of projection.
%
%  lat        - Actual lat value following 'lat'.
%
% 'trans'     - Alters transparency of the globe.
%
%  trans      - Follows 'trans' with transparency properties.
%               If set to 2, this places an opaque layer deep 
%               inside the Earth, and 3 sets the layer even deeper,
%               if set to 4, the Earth is completely transparent
%	        {Default is 1}.
%
% 'noland'    - Stop land mask being plotted.
%
% 'grid'      - Overlays a grid on the map.
% 
%
%
% OUTPUT:
%
% maph - this function generates a handle for the map
%
% A Spherical Projection is used.  Typing 'help globe' 
% in matlab will provide more information about this 
% projection.   Use the function ridgepack_sunlightm to add
% sunlight to the globe for a given date and time.
%
% Ridgepack Version 1.0
% Andrew Roberts, Naval Postgraduate School, March 2018 (afrobert@nps.edu)


global debug;
if debug; disp(['Entering ',mfilename,'...']); end

% defaults
vlat=90;
vlon=0;
land=1;
transparent=1;
grid=0;

% obtain input arguments
if nargin == 0
 disp('Default Globe');
else
 i=0;
 while i < nargin
  i=i+1;
  switch lower(varargin{i})
   case 'lat'
           i=i+1;
           if(~isnumeric(varargin{i}) || length(varargin{i})>1)
             error('Latitude must be one number')
           end
           vlat=varargin{i};
   case 'lon'
           i=i+1;
           if(~isnumeric(varargin{i}) || length(varargin{i})>1)
             error('Longitude must be one number')
           end
           vlon=varargin{i};
   case 'trans'
           i=i+1;
           if(~isnumeric(varargin{i}) || varargin{i}<1 || varargin{i}>4)
             error('Trans has the wrong value')
           end
           transparent=varargin{i};
   case 'noland'
           land=0;
   case 'grid'
           grid=1;
   otherwise
           error(['Option ',varargin{i},' is incorrect'])
  end
 end
end

% edge, grid and label color
gridcolor=0.4*[1 1 1] ; 

% clear axes
cla; 

% reset to map defaults
defaultm; 

% map projection
maph=axesm('MapProjection','globe','Frame','off');

axis off
axis vis3d

% map grid
if grid
 setm(maph,'Grid','on',...
  'GColor',gridcolor,...
  'GLinestyle',':',...
  'Glinewidth',0.5,...
  'MLinefill',5000,...
  'MLineLocation',[30],...
  'MLineLimit',[-60 60],...
  'MLineException',[-90 0 90 180],...
  'MLineVisible','on',...
  'PLineException',[],...
  'PLineFill',5000,...
  'PLineLimit',[-180 180],...
  'PLineLocation',[30],...
  'PLineVisible','on');
end

% Choose viewing angle
view(vlon+90,vlat)

% make sphere on opaque light grey, setting the rendered light color
% below the map and data, to avoid overwrite problems
if transparent==1
 base = -0.01*ones(180,361); baseref = [1 90 0];
 hs = meshm(base,baseref,size(base),base);
 set(findobj(gca,'type','surface'),'edgecolor','none','facecolor',[.97 .97 .97]);
elseif transparent==2
 base = -0.1*ones(180,361); baseref = [1 90 0];
 hs = meshm(base,baseref,size(base),base);
 set(findobj(gca,'type','surface'),'edgecolor','none','facecolor',[.97 .97 .97]);
elseif transparent==3
 base = -1*ones(180,361); baseref = [1 90 0];
 hs = meshm(base,baseref,size(base),base);
 set(findobj(gca,'type','surface'),'edgecolor','none','facecolor',[.97 .97 .97]);
else
 disp('Transparent globe selected')
end

% highlight land
if land; 
 landm = shaperead('landareas.shp', 'UseGeoCoords', true);

 % land outline (positioned at a slight altitude to remain visible)
 linem([landm.Lat],[landm.Lon],0.00001,'Color',0.45*[1 1 1],'LineWidth',0.25);
end

% remove surrounding axis
axis off
axis vis3d

% remove output if not requested
if nargout==0; clear maph; end

% Remove meridians and parallels from legend (if one is used)
h1=handlem('parallel');
hCGroup=hggroup;
set(h1,'Parent',hCGroup)
set(get(get(hCGroup,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');

h2=handlem('meridian');
hSGroup=hggroup;
set(h2,'Parent',hSGroup)
set(get(get(hSGroup,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');

drawnow

if debug; disp(['...Leaving ',mfilename]); end

