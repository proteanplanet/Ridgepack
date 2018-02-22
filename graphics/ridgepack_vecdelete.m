function ridgepack_vecdelete(ha)

% ridgepack_vecdelete - Delete  the vector key when axes change position
%
% function ridgepack_vecdelete(ha)
%
% This function deletes a vector key from a plot after it has been plotted
% with ridgepack_quiver or ridgepack_quiverm. This is useful if one wants to plot multiple
% axes on a single page, but only show a single vector key that applies to
% all panels (axes) of the figure window.
%
% Inputs:
% ha - axis handle for the quiver plot
%
% Andrew Roberts, Naval Postgraduate School, March 2018  (afrobert@nps.edu)
%

global debug;
if debug; disp(['Entering ',mfilename,'...']); end

% sub in gca if needed
if nargin==1 & ~ishandle(ha)
 error('ha is not a handle')
elseif nargin==0
 ha=gca;
end

% Extract source information
try
 h=getappdata(ha);
catch
 error('getting vector app data from axis handle')
end

% Find axis 
if isfield(h,'ReferenceVectorDeleted') & h.ReferenceVectorDeleted
 disp('Reference vector has already been deleted')
 return
elseif ~isfield(h,'ReferenceVectorStructure') | ~isfield(h,'ReferenceVectorListener') 
 error('Cannot find reference vector information')
end

% Find reference vector data structure in axis handle
rc=h.ReferenceVectorStructure;
hdl=h.ReferenceVectorListener;
dcc=h.ReferenceVectorDeleted;

% Remove text, patch, line and listener to delete reference vector
if isempty(rc) | isempty(hdl) | dcc
 error('Reference Vector Key: No Structure or Listener found')
else
 delete(rc.ht);
 delete(rc.hp);
 delete(rc.hl);
 delete(hdl); 
 rmappdata(ha,'ReferenceVectorStructure')
 rmappdata(ha,'ReferenceVectorListener')
 setappdata(ha,'ReferenceVectorDeleted',true);
end

if debug; disp(['...Leaving ',mfilename]); end


