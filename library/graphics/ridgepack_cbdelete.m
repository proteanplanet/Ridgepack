function ridgepack_cbdelete(h)

% ridgepack_cbdelete - Delete the colorbar associated with main axes
%
% function ridgepack_cbdelete(h)
%
% This function deletes an nc colorbar associated with axes handle 
% h for a figure.  Note that the axes handle h is not the handle
% of the colorbar itself, but rather of the axes whose colors
% it describes. This is analagous to using "colorbar('off')" for
% a standard matlab colorbar. 
% 
% INPUT:
%
% h - axes handle for main axes with which the colorbar is associated.
%     This can be omitted if the current axis handle is the same has h.
%
%
% OUTPUT:
%
% Output is only graphical for this function, with changes to the 
% axes handle h going on under the hood.
%
% Ridgepack Version 1.0
% Andrew Roberts, Naval Postgraduate School, March 2018 (afrobert@nps.edu)

global debug;
if debug; disp(['Entering ',mfilename,'...']); end

% check for input applicability
if nargin <1 
 h=gca;
elseif ~ishandle(h)
 error('Input is incorrect')
end

% Extract colorbar information from axis app data
hcb=getappdata(h,'ColorbarHandle');
hdl=getappdata(h,'ColorbarListener');
orientation=getappdata(h,'ColorbarOrientation');

% Delete colorbar axis handle
if ishandle(hcb)
 if ~strcmp(get(hcb,'Tag'),'Colorbar')
   error('First handle is not for an icepack colorbar')
 elseif strcmp(orientation,'vertical')
   orientation='vertical';
   rmappdata(h,'ColorbarOrientation')
   cbpos=get(hcb,'Position');
   delete(hcb)
 elseif strcmp(orientation,'horizontal')
   orientation='horizontal';
   rmappdata(h,'ColorbarOrientation')
   cbpos=get(hcb,'Position');
   delete(hcb)
 else
   delete(hcb)
   error('No colorbar orientation information found')
 end
else
 cbpos=zeros([4 1]);
end


% Delete listener for resizing colorbar
if isfield(h,'ColorbarListener');
 delete(hdl)
 rmappdata(h,'ColorbarListener')
end


% Reposition plot into center of frame, using dimensions set
% in private/ridgepack_cbpos.m to remove the colorbar width as well
% as the separation between the main axes and colorbar in pixels.

axpos=get(h,'Position');
if strcmp(orientation,'vertical')
   axpos(1)=axpos(1)+cbpos(3)/2;
elseif strcmp(orientation,'horizontal')
   axpos(2)=axpos(2)-cbpos(2)/2;
end
set(h,'Position',axpos,'Box','on');
 
axpos=getpixelposition(h);
if strcmp(orientation,'vertical')
   axpos(1)=axpos(1)+axpos(3)*0.05/2;
elseif strcmp(orientation,'horizontal')
   axpos(2)=axpos(2)+axpos(4)*0.05/2;
end
setpixelposition(h,axpos);


% Tag the app data
setappdata(h,'ColorbarHandle','ColorbarDeleted')

drawnow;

if debug; disp(['...Leaving ',mfilename]); end

