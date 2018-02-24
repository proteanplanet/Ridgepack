function ridgepack_multicb(h)

% ridgepack_multicb - Migrate a colorbar created with ridgepack_colorbar from axes to multiplot
%
% function ridgepack_multicb(h)
%
% This function moves a colorbar associated with an axes and aligns it with 
% multiple axes plotted on a single figure window.  In order to use this 
% function, you must first create axes using the ridgepack_multiplot function. 
% The end result is that instead of having a small colorbar aligned with an
% individual axes on a figure, the colorbar is moved so that it sits
% down the entire right hand side of multiple axes, in the case of a 
% vertical colorbar, or along the bottom of the figure window in the case
% of a horizontal colorbar. 
%
% Please note that the colorbar must have been created using the ridgepack_colorbar 
% function within icepack. ridgepack_multicb will not work with a plain Matlab colorbar.
% 
% INPUT:
%
% h - handle for the axes with which a colorbar is associated.
% 
% Example:
% This section provides a complete example of how ridgepack_multicb can be used 
% with icepack functions to build a 3x3 axes figure with a colorbar extending
% down the left hand side of the entire figure. The script provided here
% uses a netcdf file named "r30RB1FIY.cpl.ha.a2xavg_Sa_pslv.1992.nc"
% with x and y axes 'a2xavg_nx' and 'a2xavg_ny' to provide a color image of 
% the field 'a2xavg_Sa_pslv', respectively.  In reality, one would wish
% to use a different set of data for each figure, so that new data would need
% to be loaded for each call to ridgepack_image, and contour limits would need to be 
% provided to ridgepack_image so that each figure has a uniform color scale.
%
% clear; clf
% nc=ridgepack_clone('r30RB1FIY.cpl.ha.a2xavg_Sa_pslv.1992.nc',{'a2xavg_Sa_pslv'},1)
% ridgepack_multiplot(3,3,1,1,'a'), ridgepack_image(nc,'a2xavg_nx','a2xavg_ny','a2xavg_Sa_pslv');
% ridgepack_clearax('x',1); title(''); ridgepack_cbdelete(gca)
% ridgepack_multiplot(3,3,1,2,'b'), ridgepack_image(nc,'a2xavg_nx','a2xavg_ny','a2xavg_Sa_pslv');
% ridgepack_clearax('both',1); title(''); ridgepack_cbdelete(gca)
% ridgepack_multiplot(3,3,1,3,'c'), ridgepack_image(nc,'a2xavg_nx','a2xavg_ny','a2xavg_Sa_pslv');
% ridgepack_clearax('both',1); title(''); ridgepack_cbdelete(gca)
% ridgepack_multiplot(3,3,2,1,'d'), ridgepack_image(nc,'a2xavg_nx','a2xavg_ny','a2xavg_Sa_pslv');
% ridgepack_clearax('x',1); title(''); ridgepack_cbdelete(gca)
% ridgepack_multiplot(3,3,2,2,'e'), ridgepack_image(nc,'a2xavg_nx','a2xavg_ny','a2xavg_Sa_pslv');
% ridgepack_clearax('both',1); title(''); ridgepack_cbdelete(gca)
% ridgepack_multiplot(3,3,2,3,'f'), ridgepack_image(nc,'a2xavg_nx','a2xavg_ny','a2xavg_Sa_pslv');
% ridgepack_clearax('both',1); title(''); ridgepack_cbdelete(gca)
% ridgepack_multiplot(3,3,3,1,'g'), ridgepack_image(nc,'a2xavg_nx','a2xavg_ny','a2xavg_Sa_pslv');
% title(''); ridgepack_cbdelete(gca)
% ridgepack_multiplot(3,3,3,2,'h'), ridgepack_image(nc,'a2xavg_nx','a2xavg_ny','a2xavg_Sa_pslv');
% ridgepack_clearax('y',1); title(''); ridgepack_cbdelete(gca)
% ridgepack_multiplot(3,3,3,3,'i'), ridgepack_image(nc,'a2xavg_nx','a2xavg_ny','a2xavg_Sa_pslv');
% ridgepack_clearax('y',1); title(''); 
% ridgepack_multicb(gca)
% ridgepack_multialign(gcf,'This is the title')
%
% In reality, one would wish to use a different set of data for each figure, 
% so that new data would need to be loaded for each call to ridgepack_image, and contour 
% limits would need to be provided to ridgepack_image so that each figure has a uniform 
% color scale. This has not been done in this example for brevity.
% 
% Ridgepack Version 1.0
% Andrew Roberts, Naval Postgraduate School, March 2018 (afrobert@nps.edu)


global debug;
if debug; disp(['Entering ',mfilename,'...']); end

if nargin<1
 h=gca;
elseif ~ishandle(h)
 error('Input is not a handle')
elseif ~strcmpi(get(h,'type'),'axes')
 error('handle provided is not an axes handle')
end

% Set current figure handle
hf=gcf;

% get multifigure global axes handle
r=getappdata(hf,'MultiplotConfiguration');
if isempty(r)
 error('No multiplot configuration set for the current axes')
elseif ~isfield(r,'multihandle')
 error('handle of multiplot global axes not found')
elseif ~ishandle(r.multihandle)
 error('multiplot global axes handle error')
end

% align all before moving further
ridgepack_multialign(hf)

% get colorbar information for current axes
hcb=getappdata(h,'ColorbarHandle');
hdl=getappdata(h,'ColorbarListener');
orientation=getappdata(h,'ColorbarOrientation');

if ~isempty(hcb) & ~isempty(hdl) & ~isempty(orientation)

 % replace listener
 delete(hdl)
 hdl=addlistener(r.multihandle,{'XLim','YLim','ZLim','Position',...
                   'DataAspectRatio','PlotBoxAspectRatio',...
                   'DataAspectRatioMode','PlotBoxAspectRatioMode',...
                   'XTickLabel','YTickLabel','ZTickLabel',...
                   'OuterPosition'},'PostSet',@ridgepack_cbfix);

 % reassign to the global axes
 setappdata(r.multihandle,'ColorbarHandle',hcb)
 setappdata(r.multihandle,'ColorbarListener',hdl)
 setappdata(r.multihandle,'ColorbarOrientation',orientation)

 % accommodate colorbar
 set(r.multihandle,'ActivePositionProperty','Position')
 set(r.multihandle,'Tag','MainAxes')

 % transfer details
 rmappdata(h,'ColorbarHandle')
 rmappdata(h,'ColorbarListener')
 rmappdata(h,'ColorbarOrientation')
 set(h,'DeleteFcn','')

 % reposition colorbar
 ridgepack_cbpos(r.multihandle,orientation);

 % realign axes
 ridgepack_multialign(hf)

else

 error('No colorbar found for this axis')

end

drawnow

if debug; disp(['...Leaving ',mfilename]); end

