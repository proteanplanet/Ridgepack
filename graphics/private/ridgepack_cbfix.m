function ridgepack_cbfix(obj,src)

% function ridgepack_cbfix(obj,src)
%
% This function is part of Ridgepack Version 1.0.
% It is a listener call back to reposition the location of colorbars
% whenever there is a change to the axes accompanying the colorbar.  Note that
% this is for colorbars generated with ridgepack_colorbar, not with the vanilla colorbar
% routine available in matlab. 
% 
% Inputs:
% obj - object being called
% src - information about the affected object, which includes the 
%       affected axes handle that has been repositioned.
%
% Andrew Roberts, Naval Postgraduate School, March 2018 (afrobert@nps.edu)

% Find axis 
if ~(ishandle(src.AffectedObject) && ...
   any(strcmp(get(src.AffectedObject,'Type'),{'axes'})))
 error('src information does not contain axes handle')
end

% Call repositioning function
ridgepack_cbpos(src.AffectedObject);

