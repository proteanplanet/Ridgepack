function [NewPosition]=ridgepack_axpos(h)

% function [NewPosition]=ridgepack_axpos(h)
%
% This function is part of Ridgepack Version 1.0.
% It calculates the 'true' axes position, as seen on a figure, 
% as against the position returned by get(gca,'Position') for current axes. 
% The two may differ if the aspect ration of the data or plot box are not 1:1,
% or if the position ratio is also not square. The result is that if the
% get(gca,'Position') is used to align colorbars and plots, then large gaps
% can appear between the colorbar and its associated axes, or between axes
% on a multiplot.  As a result, this function is used when generating
% colorbars in ridgepack using ridgepack_colorbar.
%
% INPUT:
%
% h - handle of the axes for which the true edge position is required.
%
%
% OUTPUT:
%
% NewPosition - 'True' visual position of axes in normalized units
%
% Andrew Roberts, Naval Postgraduate School, March 2018 (afrobert@nps.edu)

if nargin==0
 error('No axis handle input')
elseif ishandle(h) && any(strcmp(get(h,'Type'),{'axes'}))
 axhandle=get(h);
else
 error('Input is not an axes handle')
end

% Find Position in pixels of main axes and get the tight inset
Position=getpixelposition(h);

% Caculate position ratio of main axis
posratio=Position(4)/Position(3);

% Find mode of operation setting plot limits
if strcmp(axhandle.DataAspectRatioMode,'manual') && strcmp(axhandle.PlotBoxAspectRatioMode,'manual')
 figratio=(min(diff(axhandle.YLim),axhandle.PlotBoxAspectRatio(2))*axhandle.DataAspectRatio(1))/...
          (min(diff(axhandle.XLim),axhandle.PlotBoxAspectRatio(1))*axhandle.DataAspectRatio(2));
elseif strcmp(axhandle.DataAspectRatioMode,'manual')
 figratio=(diff(axhandle.YLim)*axhandle.DataAspectRatio(1))/...
          (diff(axhandle.XLim)*axhandle.DataAspectRatio(2));
elseif strcmp(axhandle.PlotBoxAspectRatioMode,'manual')
 figratio=axhandle.PlotBoxAspectRatio(2)/...
          axhandle.PlotBoxAspectRatio(1);
else
 figratio=NaN;
end

% Find current 'real' position of axes based on position data
if isnan(figratio)

 NewPosition(3)=Position(3);
 NewPosition(1)=Position(1);

 NewPosition(4)=Position(4);
 NewPosition(2)=Position(2);

elseif posratio>=figratio

 NewPosition(3)=Position(3);
 NewPosition(1)=Position(1);

 NewPosition(4)=Position(3)*figratio;
 NewPosition(2)=Position(2)+(Position(4)-NewPosition(4))/2;

elseif posratio<figratio

 NewPosition(3)=Position(4)/figratio;
 NewPosition(1)=Position(1)+(Position(3)-NewPosition(3))/2;

 NewPosition(4)=Position(4);
 NewPosition(2)=Position(2);

end

