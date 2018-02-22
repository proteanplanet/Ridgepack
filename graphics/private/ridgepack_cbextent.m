function [tightinset]=ridgepack_cbextent(ha)

% ridgepack_cbextent - Calculate the tightinset margins of colorbars created with ridgepack_colorbar
%
% function [tightinset]=ridgepack_cbextent(ha)
%
% This function is used internally by icepack to calculate the text extent, in 
% normalized coordinates, of the colorbar scale and unit indicator.  The text
% extent is first assigned to a structure for each text marker on the colorbar. 
% This is done in ridgepack_colorbar, and is written in local data units on the axis. 
% This function then takes these extents for each entry, which can be 'scale',
% for number written on a colorbar, 'units' for the units text, or 'arrow' for
% arrows indicating colorbar continuation, and converts them onto normalized
% units on the figure window using the pixel extent and x and y limits of the 
% colorbar local data units. Please note that the arrow extents are not used
% in this function, and this function is only for use with colorbars created
% with the function "ridgepack_colorbar.m" in icepack. 
%
% Input:
% ha - handle of the main axes with which the colorbar is associated
%
% Output:
% TightInset - The tight inset of the colorbar axes, including markers, which
%              is equivalent to the TightInset obtained from normal axes. 
%
% Note:
% The original version of this function called up the extent in local data
% units, but this caused unpredictable results due to internal MATLAB
% computations and pixel rounding.  Therefore icepack was changed to immediately
% obtain the text extent when a colorbar is created, and to convert these
% to normalized coordinates whenever the colorbar position is changed.
% This resulted in a stable implementation.
% 
% Andrew Roberts, Naval Postgraduate School, March 2018  (afrobert@nps.edu)
%

global debug;
if debug; disp(['Entering ',mfilename,'...']); end

% Set the extra amount to be added to text extent for reading clarity
% This may also be required if there are negative numbers in the colorbar
% to avoid right justified colorbars from cutting off the scale.
hoverflow=1.5;
voverflow=1.5;

if ~isempty(ha) & ishandle(ha)

 hcb=getappdata(ha,'ColorbarHandle');
 orientation=getappdata(ha,'ColorbarOrientation');
 hf=get(ha,'Parent');
 figpos=getpixelposition(hf);

 if ~isempty(hcb) & ishandle(hcb) & ~isempty(orientation)
  
  tightinset=zeros([4 1]);

  axpos=getpixelposition(hcb);
  xlim=get(hcb,'XLim');
  ylim=get(hcb,'YLim');
  xextent=abs(diff(xlim));
  yextent=abs(diff(ylim));

  p=get(hcb,'Children');

  if ~isempty(p)

   widths=zeros(size(p));
   heights=zeros(size(p));
   widthu=zeros(size(p));
   heightu=zeros(size(p));
 
   for i=1:length(p)
    typep=get(p(i),'type');
    userp=get(p(i),'UserData');
    if isstruct(userp) & isfield(userp,'type') & isfield(userp,'extent')
     if strcmpi(typep,'text') & strcmpi(userp.type,'scale') 
       widths(i)=hoverflow*(userp.extent(3)/xextent)*(axpos(3)/figpos(3));
       heights(i)=voverflow*(userp.extent(4)/yextent)*(axpos(4)/figpos(4));
     elseif strcmpi(typep,'text') & strcmpi(userp.type,'units')
       if strcmp(orientation,'vertical') 
        widthu(i)=(hoverflow*(userp.extent(3)/xextent)-1)*(axpos(3)/figpos(3));
       else
        widthu(i)=hoverflow*(userp.extent(3)/xextent)*(axpos(3)/figpos(3));
       end
       heightu(i)=voverflow*(userp.extent(4)/yextent)*(axpos(4)/figpos(4));
     end
    elseif isstruct(userp)
     error('Missing fields ''type'' or ''extent'' in colorbar text descriptors')
    end
   end

   if strcmp(orientation,'vertical')
    tightinset(2)=tightinset(2)+max(heightu(:));
    tightinset(3)=tightinset(3)+max([widths(:); widthu(:)]);
    tightinset(4)=tightinset(4)+max(heightu(:))/2;
   elseif strcmp(orientation,'horizontal')
    tightinset(1)=tightinset(1)+max(widthu(:))/2;
    tightinset(2)=tightinset(2)+max(heights(:));
    tightinset(3)=tightinset(3)+max(widthu(:));
   else
    error('Colorbar orientation specification incorrect')
   end

  end

  set(hcb,'UserData',tightinset);

 elseif ~isempty(hcb) & ishandle(hcb) & isempty(orientation)
  error('Accessing colorbar orientation')
 else
  tightinset=[];
 end
 
else
  error('Problem accessing axes handle')
end

if nargout==0; clear tightinset; end

if debug; disp(['...Leaving ',mfilename]); end

