function ridgepack_trackm(nc,color,label)

% ridgepack_trackm - Draws buoy tracks on a map
%
% function ridgepack_trackm(nc,color,label)
%
% This function draws buoy tracks on a map with a label on the upper
% left side of the track, and in the color provided. A circle is plotted
% at the end of the track provided.
%
% INPUT:
%
% nc    - nc structure containing the fields latitude and longitude of the 
%         buoy track timeseries.
% color - matlab color, e.g. 'k' for black, 'b' for blue and so on
% label - text label to be placed on the upper right of the track.
%
% Ridgepack Version 1.0
% Andrew Roberts, Naval Postgraduate School, March 2018 (afrobert@nps.edu)
%

global debug;
if debug; disp(['Entering ',mfilename,'...']); end

if ~ismap(gca)
 error('Must be applied to a current map handle')
end

plotm(nc.latitude.data,nc.longitude.data,color);
plotm(nc.latitude.data(end),nc.longitude.data(end),[color,'o']);
[xpos,ypos]=mfwdtran(nc.latitude.data,nc.longitude.data);
xlab=max(xpos);
ylab=max(ypos);
text(xlab,ylab,label,'Color',color,'VerticalAlignment','bottom','HorizontalAlignment','left');

drawnow

if debug; disp(['...Leaving ',mfilename]); end


