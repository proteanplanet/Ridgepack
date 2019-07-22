function [nc]=ridgepack_hydrologym(colorindex)

% ridgepack_hydrologym - Colors location of lakes on a map.
%
% function ridgepack_hydrologym(colorindex)
%
% Colors location of lakes blue for the current axes handle.
% Also adds in prominant rivers to the current map. There is a
% choice of four different colors for the rivers.
%
% INPUT:
%
% colorindex - integer from 1 to 4 that selects the desired 
%              color (optional)
%
% Ridgepack Version 1.0
% Andrew Roberts, Naval Postgraduate School, March 2018 (afrobert@nps.edu)
 

global debug;
if debug; disp(['Entering ',mfilename,'...']); end

if ~ismap(gca)
 error('Must be applied to a current map handle')
end

if nargin>=1
 if colorindex==1
  truecolor=0.5*[0.8 1 0.8];
 elseif colorindex==2
  truecolor=0.75*[0 1 0];
 elseif colorindex==3
  truecolor=[0.75 1 1];
 elseif colorindex==4
  truecolor=0.9*[0 0 1];
 else
  error('Color index chosen is out of range')
 end
else
 truecolor=0.9*[0 0 1];
end

lakes = shaperead('worldlakes','UseGeoCoords',true);
geoshow(gca,lakes, 'FaceColor',truecolor,'EdgeColor',truecolor,'FaceAlpha',0.25)

rivers = shaperead('worldrivers','UseGeoCoords',true);
geoshow(gca,rivers, 'Color',truecolor,'LineWidth',0.8,'DisplayType','line')

nc.attributes.title='Global River Tracks';
nrivers=length(rivers);
maxlength=0;
for i=1:nrivers
 maxlength=max(maxlength,length(rivers(i).Lat));
end

nc.nrivers.data=[1:nrivers];
nc.nrivers.long_name='River Number Index';
nc.nrivers.dimension={'nrivers'};
nr.nrivers.units='';

nc.points.data=[1:maxlength];
nc.points.long_name='Point on trace';
nc.points.dimension={'points'};
nr.points.units='';

nc.latitude.data=NaN*zeros([nrivers maxlength]);
nc.longitude.data=NaN*zeros([nrivers maxlength]);
for i=1:nrivers
 riverlength=length(rivers(i).Lat);
 nc.latitude.data(i,1:riverlength)=rivers(i).Lat;
 nc.longitude.data(i,1:riverlength)=rivers(i).Lon;
end

nc.latitude.long_name='Latitude of Rivers';
nc.latitude.dimension={'nrivers','points'};
nc.latitude.units='degrees_north';

nc.longitude.long_name='Longitude of Rivers';
nc.longitude.dimension={'nrivers','points'};
nc.longitude.units='degrees_east';

nc=ridgepack_struct(nc);

drawnow

if debug; disp(['...Leaving ',mfilename]); end

