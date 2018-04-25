function ridgepack_hydrologym(colorindex)

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
geoshow(gca,rivers, 'Color',truecolor,'LineWidth',0.4,'DisplayType','line')

drawnow

if debug; disp(['...Leaving ',mfilename]); end

