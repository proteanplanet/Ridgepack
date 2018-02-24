function ridgepack_cbfix(obj,src)

% ridgepack_cbfix - Callback to reposition colorbar when axes change position
%
% function ridgepack_cbfix(obj,src)
%
% This function is a listener call back to reposition the location of colorbars
% whenever there is a change to the axes accompanying the colorbar.  Note that
% this is for colorbars generated with ridgepack_colorbar, not with the vanilla colorbar
% routine available in matlab. 
% 
% INPUT:
%
% obj - object being called
% src - information about the affected object, which includes the 
%       affected axes handle that has been repositioned.
%
% Ridgepack Version 1.0
% Andrew Roberts, Naval Postgraduate School, March 2018 (afrobert@nps.edu)

global debug;

if debug; disp(['Entering ',mfilename,'...']); end

% Find axis 
if ~(ishandle(src.AffectedObject) && ...
   any(strcmp(get(src.AffectedObject,'Type'),{'axes'})))
 error('src information does not contain axes handle')
end

% Call repositioning function
ridgepack_cbpos(src.AffectedObject);

if debug; disp(['...Leaving ',mfilename]); end

