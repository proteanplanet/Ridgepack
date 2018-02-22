function ridgepack_cbshare(h)

% ridgepack_cbshare - Share a colorbar from one axes with adjacent axes
%
% function ridgepack_cbshare(h)
%
% This function shares a vertical colorbar from one axes with the 
% axes immediately beneath it, designed for when multiple axes are 
% plotted in a single figure.  For horizontal axes, the colorbar
% is shared with the axes immediately to the right of the axes in 
% a multiple axes figure. What this means is that the colorbar
% is either shifted down or to the right by half the axes width
% so that visually the colorbar spans a set of two vertical or
% horizontal axes.  This can be used when a colorbar in a multiplot
% figure only applies to two rows or two columns of axes, rather 
% than for every frame.  If the colorbar applies to all axes, then
% ridgepack_multicb should be used instead.  This function only works
% for colorbars created with ridgepack_colorbar, not with generic MATLAB
% colorbars.
%
% Input:
% h - main axes with the colorbar to be shared.
%
% Andrew Roberts, Naval Postgraduate School, March 2018  (afrobert@nps.edu)
%

% Find axis and colorbar handles
if nargin==0
 error('No axis handle input')
elseif ishandle(h) && any(strcmp(get(h,'Type'),{'axes'}))
 axhandle=get(h);
 hcb=getappdata(h,'ColorbarHandle');
else
 error('Input is not an axes handle')
end

% Set sharing to be turned on
if ishandle(hcb)
 setappdata(h,'ColorbarShare',1);
 orientation=getappdata(h,'ColorbarOrientation');
 if strcmp(orientation,'vertical')
  disp('Sharing colorbar with axes beneath this one')
 elseif strcmp(orientation,'horizontal')
  disp('Sharing colorbar with axes to the right of this one')
 else
  error('Incorrect orientiation specification')
 end
elseif ischar(hcb) && strcmp(hcb,'ColorbarDeleted')
 % reached here if colorbar has been deleted 
 return
else
 error('Unable to share the colorbar')
end

drawnow


