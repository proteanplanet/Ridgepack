function [mstruct]=ridgepack_extendframe

% ridgepack_extendframe - This function extends the frame around map to limit trimming
%
% function [mstruct]=ridgepack_extendframe
%
% This program should be run when there is a map open and active.  It uses
% the current map structure and and extends the bounds of the frame in a 
% new identical mapstructure, mstruct, which can then be used for translating
% data to plotting coordinates from geographical coordinates.  
%
% This function is necessary because for many map projections, regardless
% of the frame set in the map setup, it is contracted to the edge of the plotted
% map.  This means that data outside this region is trimmed. Therefore grid
% cells only fractionally in the frame are removed, resulting in vectors, contours
% and patches not extending to the very edge of the map being plotted.  
%
% The frame is extended past the edge, where possible, by 5 degrees, which 
% will typically take 2.5 degree datasets two grid cells past the edge of the
% frame.
%
% OUTPUT:
%
% mstruct - map structure
%
% Ridgepack Version 1.0
% Andrew Roberts, Naval Postgraduate School, March 2018 (afrobert@nps.edu)


global debug;
if debug; disp(['Entering ',mfilename,'...']); end

extend=5; % degrees

mstruct=gcm;

if strcmp(mstruct.mapprojection,'stereo') | strcmp(mstruct.mapprojection,'gstereo')

 mstruct.flatlimit(1)=max(mstruct.flatlimit(1)-extend,-90);
 mstruct.flatlimit(2)=min(mstruct.flatlimit(2)+extend,90);
 mstruct.flonlimit(1)=max(mstruct.flonlimit(1)-extend,-180);
 mstruct.flonlimit(2)=min(mstruct.flonlimit(2)+extend,180);

end

if debug; disp(['...Leaving ',mfilename]); end

