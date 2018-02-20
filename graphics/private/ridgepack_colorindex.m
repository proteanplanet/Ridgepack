function [zindex,truecol]=ridgepack_colorindex(z,cont,ref,mask)

% function [zindex,truecol]=ridgepack_colorindex(z,cont,ref,mask)
%
% This function is part of Ridgepack Version 1.0.
% It generates a color index and true color array from a 
% data array to be used in various matlab color plotting utilities 
% given the original data values and the contour values, as well
% the color reference value and mask.  
%
% INPUT:
%
% z      - data to be mapped
% cont   - vector of contour values (colormap boundaries)
% ref    - reference contour about which color is referenced
% mask   - mask to be added to true colors of same size as z
%
%
% OUTPUT:
%
% zindex  - color indexed data to be used for graphics plotting.
% truecol - true color array corresponding to zindex.
%
% Andrew Roberts, Naval Postgraduate School, March 2018 (afrobert@nps.edu)

% check for limited color bands
if length(cont)<3
 warning('Only a single colorband has been specified')
end

% sort cont in ascending numerical order
cont=sort(cont);

% turn shown array into indexed colors
zindex=ones(size(z));
for i=2:(length(cont)-1)
 zindex(z>=cont(i) & z<cont(i+1))=i;
end
zindex(z>=cont(end-1))=length(cont)-1;


% set NaNs in zindex, including where ref
% is the minimum or maximum contour value
% and values are less or greater than.
zindex(isnan(z))=NaN;
if ref<=cont(1)
 zindex(z<ref)=NaN;
elseif ref>=cont(end)
 zindex(z>ref)=NaN;
end


% set zindex where there is a mask and no data
if nargin==4 && all(size(mask)==size(zindex))
 zindex(mask==0 & isnan(z))=-1;
end


% find true color array
xsize=size(zindex);
if length(xsize)>2
 disp('Unable to allocate true color array')
 truecol=NaN;
else
 xsize(3)=3;
 cmap=colormap;
 cmap(end+1,:)=[1 1 1]; 
 zindex=reshape(zindex,size(zindex(:)));
 truecol=zeros([length(zindex) 3]);
 for i=1:length(zindex)
  if isnan(zindex(i))
   truecol(i,:)=NaN; % NaNs in data are blank
  elseif zindex(i)<0
   truecol(i,:)=0.96*[1 1 1]; % masked areas are grey
   zindex(i)=NaN;
  else
   truecol(i,:)=cmap(zindex(i),:); % Allocate cmap colors
  end
 end
 zindex=reshape(zindex,xsize(1:2));
 truecol=reshape(truecol,xsize(1:3));
end

