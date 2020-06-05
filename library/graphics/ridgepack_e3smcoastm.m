function [nc,SCP,SCL]=ridgepack_e3smcoastm(ncvert,shape);

% ridgepack_e3smcoastm - draw or generate an E3SM coastline
%
% function [nc,SCP,SCL]=ridgepack_e3smcoastm(ncvert,shape) 
%
% This function generates a coast given an MPAS mesh in an
% nc structure.  If output arguments are provided, then the 
% function simply generates an nc structure as well as line and 
% polygonal geoshapes. If not, the coast is plotted on a map.
%
% INPUT:
%
% ncvert - nc structure containing the following fields from 
%          an E3SM Ocean initial condition file. It can also 
%          be the geoshape SCP generated from a previous call
%          if plotting only, so as to save time.
% shape  - logical determining if it is a shape or line to be 
%          be plotted for the coast.
%
% OUTPUT:
%
% nc  - nc structure with close-loop lines describing the coastline.
% SCP - Polygons of the E3SM coastline as a geoshape.
% SCL - Line of the E3SM coastline as a geoshape.
%
% Ridgepack Version 2.0
% Andrew Roberts, Los Alamos National Laboratory, 2020

global debug;
if debug; disp(['Entering ',mfilename,'...']); end

if nargin<1
 error('missing input')
elseif ~isstruct(ncvert)
 disp('ncvert is not an nc structure, assuming a geoshape')
 if nargout==0
  SCP=ncvert;
 else
  error('must supply a netcdf structure to generate the coast')
 end
end

if nargin<2 
 shape=true;
elseif ~islogical(shape)
 error('shape should be true or false')
end

% generate the coast
if isstruct(ncvert)
 [nc,SCP,SCL]=ridgepack_e3smseasaw(ncvert);
end

% plot the coast if required
if nargout==0

 ht=get(gcf,'CurrentAxes');

 if ~ismap(ht)
  error('Current axes must be a map')
 end

 if shape
  geoshow(SCP,'FaceColor',0.95*[1 1 1],'EdgeColor',0.25*[0 0 1]);
 else
  [c,d] = mfwdtran(gcm,[SCP.Latitude],[SCP.Longitude]);
  plot(c,d,'Color',0.25*[0 0 1],'Linewidth',0.5);
 end

end

% debug stuff
if debug; disp(['...Leaving ',mfilename]); end

