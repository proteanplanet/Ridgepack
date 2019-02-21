function [maph]=ridgepack_worldm(varargin)

% ridgepack_worldm - Generates a global miller projection with land mask
%
% function [maph]=ridgepack_worldm(varargin)
%
% This generates a miller projection of the globe. 
%
% INPUT:
%
% center - Central meridian of the plot east from
%          the Greenwich Median (between 0 and 360)
%
% noland - switch off land
%
% This function generates a handle for the map:  maph
%
% A Miller Cylindrical Projection is used.  Typing 
% 'help miller' in matlab will provide more information 
% about this projection.  To see the distortion provided 
% by this map, use the 'nctissotm' function.
%
% Example:  ncworld('center',70)
%
% Ridgepack Version 1.0
% Andrew Roberts, Naval Postgraduate School, March 2018 (afrobert@nps.edu)


global debug;
if debug; disp(['Entering ',mfilename,'...']); end

% defaults 
center=180.0;
noland=0;

% obtain input arguments
if nargin == 0
 disp('Default world map');
else
 i=0;
 while i < nargin
  i=i+1;
  switch lower(varargin{i})
   case 'center'
    i=i+1;
    if(~isnumeric(varargin{i}) || length(varargin{i})>1)
      error('Center must be one number')
    end
    center=varargin{i};
   case 'noland'
    noland=1;
   otherwise
    error(['Option ',char(varargin{i}),' is incorrect'])
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
maph=axesm('MapProjection','miller',...
  'AngleUnits','degrees',...
  'Aspect','normal',...
  'FalseNorthing',0,...
  'FalseEasting',0+center,...
  'Geoid',[1 0],...
  'MapLatLimit',[-88 88],...
  'MapLonLimit',[-180+center 180+center],...
  'Scalefactor',1,...
  'Frame','on',...
  'FFill',2000,...
  'FEdgeColor',gridcolor,...
  'FFaceColor','white',...
  'FLatLimit',[-90 90],...
  'FLonLimit',[-180 180],...
  'FLineWidth',0.25);

% map grid
%setm(maph,'Grid','on',...
%  'Galtitude',Inf,...
%  'GColor',gridcolor,...
%  'GLinestyle',':',...
%  'Glinewidth',0.5,...
%  'MLinefill',2000,...
%  'MLineLimit',[-90 90],...
%  'MLineException',[-90 0 90 180],...
%  'MLineLocation',[90],...
%  'MLineVisible','on',...
%  'PLineException',[],...
%  'PLineFill',2000,...
%  'PLineLimit',[-180+center 180+center],...
%  'PLineLocation',[30],...
%  'PLineVisible','on')

% grid labels
setm(maph,'Fontangle','normal',...
  'FontColor',gridcolor,...
  'Fontname','helvetica',...
  'FontSize',8,...
  'FontUnits','points',...
  'FontWeight','normal',...
  'LabelFormat','compass',...
  'LabelRotation','off',...
  'MeridianLabel','on',...
  'MLabelLocation',[90],...
  'MLabelParallel',-97,...
  'MLabelRound',0,...
  'ParallelLabel','on',...
  'PLabelLocation',[30],...
  'PLabelMeridian',-187+center,...
  'PLabelRound',0);

if noland==0
 % land mask outline
 landm = shaperead('landareas.shp', 'UseGeoCoords', true);

 % land mask set just below the surface
 patchm([landm.Lat],[landm.Lon],-0.00001,0.96*[1 1 1],'EdgeColor','none');

 % land outline (positioned at a slight altitude to remain visible)
 linem([landm.Lat],[landm.Lon],0.00001,'Color',0.45*[1 1 1],'LineWidth',0.25);
end

% maximum plot size
tightmap 

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

