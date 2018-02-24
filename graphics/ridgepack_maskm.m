function [h]=ridgepack_maskm(a,b,zi,co,wi,li,alt)

% ridgepack_maskM - Adds a border between 1s and 0s on a mask
%
% function [h]=ridgepack_maskm(a,b,zi,co,wi,li,alt)
% 
% This function draws a black border between regions of 1s and 0s
% on a map.  It is specifically intended to highlight the land-sea 
% boundaries in numerical model output from climate, weather and ocean 
% models or for highlighting regions of diffent numerical value.
%
% INPUT:
%
% a  - latitude coordinates
% b  - longitude coordinates
% zi - values of array containing values as a mask.  This should
%      be a 2D array of 1s and 0s.
% co - optional color (RGB matrix or specific color, default is grey)
% wi - width of line (default is 0.75)
% li - linestyle (e.g. ':', default is '-').  Note that this will not
%      work for strongly curvilinear cases.
% alt- altitude of point in km
% 
%
% OUTPUT:
%
% h  - handle for generating graphic legends
%
% Ridgepack Version 1.0
% Andrew Roberts, Naval Postgraduate School, March 2018 (afrobert@nps.edu)


global debug;
if debug; disp(['Entering ',mfilename,'...']); end

disp('Drawing mask or divider on the map');

if ~ismap(gca)
 error('Must be applied to a current map handle')
end


if nargin<4; co=[0.5 0.5 0.5]; end
if nargin<5; wi=0.75; end
if nargin<6; li='-'; end

% Draw mask around 0s,keeping in mind that the 
% x-coordinate corresponds to the j index of z(i,j), 
% and the y-coordinate corresponds to the i index
% of x(i,j).

if not(ndims(zi)==2);error('zi is not a 2D array');end
if size(a)~=size(zi);error('latitude not same size as zi');end 
if size(b)~=size(zi);error('longitude not same size as zi');end 


% DO NOT DO THIS, AS IT TURNS ON BEHAVIOR IN OTHER ROUTINES NOT 
% DESIRED...
% Include NaN Values
%if any(isnan(zi(:))) & range(zi(:))~=1
% zi(~isnan(zi))=1; 
% zi(isnan(zi))=0; 
%end

% extend the arrays to account for boundaries
z(1:size(zi,1)+2,1:size(zi,2)+2)=0;
z(2:size(zi,1)+1,2:size(zi,2)+1)=zi(:,:);

x=z;
y=z;

x(2:size(zi,1)+1,2:size(zi,2)+1)=a(:,:);
x(1,2:size(zi,2)+1)=a(1,:);
x(end,2:size(zi,2)+1)=a(end,:);
x(2:size(zi,1)+1,1)=a(:,1);
x(2:size(zi,1)+1,end)=a(:,end);
x(1,1)=a(1,1);
x(1,end)=a(1,end);
x(end,1)=a(end,1);
x(end,end)=a(end,end);

y(2:size(zi,1)+1,2:size(zi,2)+1)=b(:,:);
y(1,2:size(zi,2)+1)=b(1,:);
y(end,2:size(zi,2)+1)=b(end,:);
y(2:size(zi,1)+1,1)=b(:,1);
y(2:size(zi,1)+1,end)=b(:,end);
y(1,1)=b(1,1);
y(1,end)=b(1,end);
y(end,1)=b(end,1);
y(end,end)=b(end,end);

% get current axes and convert points to mapping cartesian coordinates
% make allowance for stupid matlab trimming rules which clip just inside 
% instead of just outside of mapping frame.
h=get(gcf,'CurrentAxes');
mapobj=gcm;

if strmatch(mapobj.mapprojection,'globe')

 disp('Plotting on 3D globe')

 % set altitude before transforming data
 alt=0*ones(size(x));
 [x,y,e] = mfwdtran(x,y,alt,h);

else

 disp('Plotting on 2D map projection')
 mstruct=ridgepack_extendframe;

 [x,y]=mfwdtran(mstruct,x,y,h,'surface'); % untrimmed values

 % set altitude after transforming data
 e=ones(size(x));

end

% now draw the dividing lines between pixels
for i=2:size(zi,1)+1;
for j=2:size(zi,2)+1;
  
  if z(i,j)==0 && z(i-1,j)==1;
    lx1=mean(reshape(x(i-1:i,j-1:j),4,1));
    ly1=mean(reshape(y(i-1:i,j-1:j),4,1));
    lz1=mean(reshape(e(i-1:i,j-1:j),4,1));
    lx2=mean(reshape(x(i-1:i,j:j+1),4,1));
    ly2=mean(reshape(y(i-1:i,j:j+1),4,1));
    lz2=mean(reshape(e(i-1:i,j:j+1),4,1));
    line([lx1 lx2],[ly1 ly2],[lz1 lz2],'LineStyle',li,'Color',co,'linewidth',wi);
  end

  if z(i,j)==0 && z(i+1,j)==1;
    lx1=mean(reshape(x(i:i+1,j-1:j),4,1));
    ly1=mean(reshape(y(i:i+1,j-1:j),4,1));
    lz1=mean(reshape(e(i:i+1,j-1:j),4,1));
    lx2=mean(reshape(x(i:i+1,j:j+1),4,1));
    ly2=mean(reshape(y(i:i+1,j:j+1),4,1));
    lz2=mean(reshape(e(i:i+1,j:j+1),4,1));
    line([lx1 lx2],[ly1 ly2],[lz1 lz2],'LineStyle',li,'Color',co,'linewidth',wi);
  end

  if z(i,j)==0 && z(i,j-1)==1;
    lx1=mean(reshape(x(i-1:i,j-1:j),4,1));
    ly1=mean(reshape(y(i-1:i,j-1:j),4,1));
    lz1=mean(reshape(e(i-1:i,j-1:j),4,1));
    lx2=mean(reshape(x(i:i+1,j-1:j),4,1));
    ly2=mean(reshape(y(i:i+1,j-1:j),4,1));
    lz2=mean(reshape(e(i:i+1,j-1:j),4,1));
    line([lx1 lx2],[ly1 ly2],[lz1 lz2],'LineStyle',li,'Color',co,'linewidth',wi);
  end

  if z(i,j)==0 && z(i,j+1)==1;
    lx1=mean(reshape(x(i-1:i,j:j+1),4,1));
    ly1=mean(reshape(y(i-1:i,j:j+1),4,1));
    lz1=mean(reshape(e(i-1:i,j:j+1),4,1));
    lx2=mean(reshape(x(i:i+1,j:j+1),4,1));
    ly2=mean(reshape(y(i:i+1,j:j+1),4,1));
    lz2=mean(reshape(e(i:i+1,j:j+1),4,1));
    h=line([lx1 lx2],[ly1 ly2],[lz1 lz2],'LineStyle',li,'Color',co,'linewidth',wi);
  end

end
end

if nargout==0; clear h; end

drawnow

if debug; disp(['...Leaving ',mfilename]); end


