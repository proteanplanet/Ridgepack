function ridgepack_conicm(varargin)

% ridgepack_conicm - Conic projections on the globe
%
% function ridgepack_conicm(varargin)
%
% This function generates a conic projection over an area of the globe.
% It is limited to preset arguments listed below.
%
% INPUT:
%
% 'noland'    - Stop land mask being plotted.
%
% 'nogrid'    - Stop grid being plotted.
%
% 'nolabel'   - Stop grid labels.
%
% 'labeloff'  - Stop printing axis labels.
%
% 'beaufort'  -	Provides a conic projection centered over the Beaufort Sea,
%		which includes a grid labelled with latitudes and longitudes.
%	        This is the default projection.
%
% 'beaufort2' -	Zoomed in section of the beaufort sea.
%
% 'bering'    -	Bering Strait.
%
% 'greenland' - Conic projection of Greenland and surrounding seas.
%
% 'labrador'  - Conic projection of the Labrador and surrounding seas.
%
% 'gulfofalaska' - Conic projection of the Gulf of Alaska.
%
% 'mertz'     - mertz glacier polynya
%
%
% OUTPUT:
%
% This function generates a handle for the map:  maph
%
% Ridgepack Version 1.0
% Andrew Roberts, Naval Postgraduate School, March 2018 (afrobert@nps.edu)


global debug;
if debug; disp(['Entering ',mfilename,'...']); end

% Check that the mapping toolbox is installed
h=ver('map') ;
if length(h)==0 ; error('Mapping toolbox not installed') ; end

% defaults
beaufort=0;
beaufort2=0;
bering=0;
greenland=0;
labrador=0;
gulfofalaska=0;
mertz=0;
centralmeridian=0 ; 
equatorextent=60 ; 
parallelstop=80;
grid=0;
label=0;
mult=1;
land=1;
nolabel=false;
nogrid=false;

% edge, grid and label color
gridcolor=[0.400 0.400 0.400] ; 

% set default map properties
defaultm; 

% option selection
if nargin == 0
 disp('Default Beaufort Sea conic projection');
 beaufort=1;
else
 i=0;
 while i < nargin
  i=i+1;
  switch lower(varargin{i})
   case 'beaufort'
	beaufort=1;
   case 'beaufort2'
	beaufort2=1;
   case 'bering'
	bering=1;
   case 'greenland'
	greenland=1;
   case 'labrador'
	labrador=1;
   case 'gulfofalaska'
	gulfofalaska=1;
   case 'mertz'
	mertz=1;
   case 'noland'
        land=0;
   case 'nolabel'
        nolabel=true;
   case 'nogrid'
        nogrid=true;
   otherwise
	error(['Option ',varargin{i},' is incorrect'])
  end
 end
end

% this is changed by only a few map configs
gridseperation=5;

if greenland
 minlat=56;
 maxlat=85;
 minlon=-75;
 maxlon=-5;
 MLabelLocation=[-70:10:-10];
 PLabelLocation=[60:5:85];
elseif labrador
 minlat=50;
 maxlat=75;
 minlon=-70;
 maxlon=-30;
 MLabelLocation=[-70:10:-30];
 PLabelLocation=[50:5:75];
elseif gulfofalaska
 minlat=28;
 maxlat=65;
 minlon=170;
 maxlon=250;
 MLabelLocation=[-360:10:360];
 PLabelLocation=[30:5:65];
elseif beaufort
 minlat=67;
 maxlat=80;
 minlon=-165;
 maxlon=-115;
 MLabelLocation=[-160:10:-120];
 PLabelLocation=[70:5:80];
elseif beaufort2
 minlat=73;
 maxlat=78;
 minlon=-150;
 maxlon=-130;
 MLabelLocation=[-160:10:-120];
 PLabelLocation=[70:5:85];
elseif bering
 minlat=60;
 maxlat=70;
 minlon=-180;
 maxlon=-155;
 MLabelLocation=[-180:5:-155];
 PLabelLocation=[60:5:70];
elseif mertz
 minlat=-68;
 maxlat=-64;
 minlon=137;
 maxlon=149;
 gridseperation=1;
 MLabelLocation=[minlon:gridseperation:maxlon];
 PLabelLocation=[minlat:gridseperation:maxlat];
else
 error('Error: Please define the area of the conic map')
end

% map projection
maph=axesm('MapProjection','eqaconic',...
  'AngleUnits','degrees',...
  'MapParallels',[minlat maxlat], ...
  'Aspect','normal',...
  'FalseNorthing',0,...
  'FalseEasting',0,...
  'Geoid',[1 0],...
  'MapLatLimit',[minlat maxlat],...
  'MapLonLimit',[minlon maxlon],...
  'FLatLimit',[minlat maxlat],...
  'FLonLimit',[minlat maxlat],...
  'Frame','off',...
  'FFill',2000,...
  'FLineWidth',eps);

hh=getm(maph);

% add land to plot
if (land==1)
 landm = shaperead('landareas.shp', 'UseGeoCoords', true);

 % land mask 
 patchm([landm.Lat],[landm.Lon],-0.1,0.96*[1 1 1],'EdgeColor','none');

 % land outline (positioned at a slight altitude to remain visible)
 linem([landm.Lat],[landm.Lon],0.00001,'Color',0.45*[1 1 1],'LineWidth',0.25);
end

% draw grid
if ~nogrid
 setm(maph,'Grid','on',...
  'Galtitude',1,...
  'GColor',gridcolor,...
  'GLinestyle',':',...
  'Glinewidth',0.5,...
  'MLinefill',3000,...
  'MLineLimit',[90 -90],... 
  'MLineException',[-90 0 90 180],...
  'MLineLocation',[gridseperation],...
  'MLineVisible','on',...
  'PLineException',[],...
  'PLineFill',2000,...
  'PLineLimit',[-360 360],...
  'PLineLocation',[gridseperation],...
  'PLineVisible','on');
end

% draw grid labels
if ~nolabel
 setm(maph,'Fontangle','normal',...
  'FontColor',[0 0 0],...
  'Fontname','helvetica',...
  'FontSize',8,...
  'FontUnits','points',...
  'FontWeight','normal',...
  'LabelFormat','none',... % could also be 'compass'
  'LabelRotation','on',...
  'MeridianLabel','on',...
  'MLabelLocation',MLabelLocation,...
  'MLabelParallel',[sign(minlat)*min(abs(maxlat),abs(minlat))],...
  'MLabelRound',0,...
  'ParallelLabel','on',...
  'PLabelLocation',PLabelLocation,...
  'PLabelMeridian',[minlon],...
  'PLabelRound',0);
end

% fit map tightly to axes
if mertz | labrador
 axis on
 framem off
else
 axis off
 framem on
end
tightmap

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

h2=handlem('Frame');
set(h2,'Clipping','on');
hSGroup=hggroup;
set(h2,'Parent',hSGroup);
set(get(get(hSGroup,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');

drawnow

limsx=xlim;
limsy=ylim;
x1=limsx(1);
x2=limsx(2);
y1=limsy(1);
y2=limsy(2);

if labrador==1
 z=0.002;
 plotbox=plot3([x1 x2 x2 x1 x1],[y1 y1 y2 y2 y1],[z z z z z],...
               'Color','k');
 uistack(plotbox,'top')
 drawnow
end

if debug; disp(['...Leaving ',mfilename]); end


