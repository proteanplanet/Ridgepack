function ridgepack_multifix(obj,src)

% ridgepack_multifix - Callback to align multiple axes on a single figure
%
% function ridgepack_multifix(obj,src)
%
% This function is a listener call back to align multiple axes on the one figure
% that have been generated with the ridgepack_multiplot function.
% 
% Inputs:
% obj - object being called
% src - information about the affected object, which includes the 
%       affected axes handle that has been repositioned.
%
% Andrew Roberts, Naval Postgraduate School, March 2018  (afrobert@nps.edu)
%

global debug;
if debug; disp(['Entering ',mfilename,'...']); end

ho=obj.Name;

if strcmp(obj.Name,'CurrentAxes')
 if isfield(src,'NewValue')
  ha=src.NewValue;
 else
  ha=1;
 end
 hf=src.AffectedObject;
else
 ha=src.AffectedObject;
 hf=ha.Parent;
end


if isempty(ha) | ~ishandle(ha)

 % remove listener and config if ha is empty or not a handle
 rmappdata(hf,'MultiplotConfiguration')
 hlf=getappdata(hf,'MultiplotFigureListener');
 delete(hlf)
 rmappdata(hf,'MultiplotFigureListener')

else

 % This commented out because it was causing a problem with the notation
 % at the top of axes

 %try
 % ridgepack_multialign(hf);
 %catch
 % % remove listener if ridgepack_multialign throws an error
 % hlf=getappdata(hf,'MultiplotFigureListener');
 % delete(hlf)
 %end

end

if debug; disp(['...Leaving ',mfilename]); end
