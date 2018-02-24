function ridgepack_vecfix(obj,src)

% ridgepack_vecfix - Callback to reposition vector key when axes change position
%
% function ridgepack_vecfix(obj,src)
%
% This function is a listener callback to reposition the vector key for
% a (non-color coded) vector plot whenever there is a change in axes position
% or limits. Whenever the axes are repositioned, the vector key is deleted,
% and a new one is redrawn in the correct place on the plot, with the new 
% handle for the key being added to the main axes app data.
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
if ~(ishandle(src.AffectedObject) && strcmp(get(src.AffectedObject,'Type'),'axes'))
 error('src information does not contain axes handle')
end

% Find reference vector data structure in axis handle
rc=getappdata(src.AffectedObject,'ReferenceVectorStructure');
hdl=getappdata(src.AffectedObject,'ReferenceVectorListener');
dcc=getappdata(src.AffectedObject,'ReferenceVectorDeleted');

% Remove text, patch, line and listener to delete reference vector
if dcc
 if debug; disp('Reference vector key: Previously Deleted'); end
 return
elseif isempty(rc) | isempty(hdl)
 if debug; disp('Reference vector key: No structure or listener'); end
else
 if debug; disp('Reference vector key: Fixing vector reference'); end
 delete(rc.ht)
 delete(rc.hp)
 delete(rc.hl)
 delete(hdl)
 rmappdata(src.AffectedObject,'ReferenceVectorStructure')
 rmappdata(src.AffectedObject,'ReferenceVectorListener')
 ridgepack_vecref(src.AffectedObject,rc.scalelength,rc.arrow,rc.refval,rc.units,...
          rc.veccol,rc.vecwidth,rc.accenture,rc.lx,rc.ly);
end

if debug; disp(['...Leaving ',mfilename]); end

