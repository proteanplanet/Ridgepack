function [x,y,u,v]=ridgepack_quiverthin(x,y,u,v,times,method)

% ridgepack_quiverTHIN - Thin vectors on a quiver plot for a map or Cartesian axis
%
% function [x,y,u,v]=ridgepack_quiverthin(x,y,u,v,times,sample)
%
% This function thins the u and v component of a quiver plot, together with their
% x and y coordinates.  The thinning applies to both coordinates in lat-long or
% in Cartesian, but note that if in lat-long, the map on which the vectors are to 
% be drawn must already be plotted ready for vector overlay, and must be the current
% axes.
%
% Input:
% x       - x-coordinate or latitude
% y       - y-coordinated or longitude
% u       - u-component (Cartesian +x-direction, map +longitude-direction)
% v       - v-component (Cartesian +y-direction, map +latitude-direction)
% times   - the number of times the thinning is to be done.  One thinning
%           reduces the vectors by a factor of two (plus removal of edge 
%           vectors where an uneven dimension size exists).  times may be
%           set to 1, to reduce the number of vectors by two, 2 to reduce 
%           the number of vectors by four, 3 to reduces the number of 
%           vectors by eight and so on.
% method  - selector indicating method for thinning vectors:
%           1: Simple sub-sampling of the field
%           2: Calculate mean, only returning values where all points are non NaNs
%           3: Calculate median with no more than half the points as NaNs to account
%              for areas next to coastlines (areas masked out with NaNs in u and v). 
%
% Output:
% x       - thinned x-coordinate or latitude
% y       - thinned y-Coordinated or longitude
% u       - thinned u-component (Cartesian +x-direction, map +longitude-direction)
% v       - thinned v-component (Cartesian +y-direction, map +latitude-direction)
%
% Andrew Roberts, Naval Postgraduate School, March 2018  (afrobert@nps.edu)
%

global debug;
if debug; disp(['Entering ',mfilename,'...']); end

% Set default values
if nargin<5 ; error('Missing input data'); end
if nargin<6; 
 method=3;
elseif ~isnumeric(method)
 error('Method (argument 6) is should be a number between 1 and 3')
elseif method<1 | method>3
 error('Method selection is correct');
end
h=get(gcf,'CurrentAxes');


% Determine if the axes are map or cartesian, if the former calculate mapping to 
% plot axis, and then thin the vector field otherwise just thin the vector field.
if ismap(h)

 disp(['Thinning vectors on current map for every ',num2str(2.^(times)),' cell(s)']);

 % set lat and lon
 lat=x;
 lon=y;

 % get x and y location on the map
 [x,y] = mfwdtran(lat,lon);

else

 disp(['Thinning vectors on cartesian axes for every ',num2str(2.^(times)),' cell(s)']);

end

if length(times)>1
 error('The vector thinning specification must be a number, not a vector')
end

mediansample=true;

if times>0 

 stepi=2.^times;
 nsamples=stepi*stepi;

 % thin vectors by taking median
 if method==3

  disp(['Thinning vectors with median using between ',num2str(nsamples/2),' and ',...
         num2str(nsamples),' samples per vector'])

  for i=1:stepi:size(x,1); for j=1:stepi:size(x,2)
   submat=x(i:min(i+stepi-1,size(x,1)),j:min(j+stepi-1,size(x,2)));
   if any(isnan(submat(:)))
    X(1+(i-1)/stepi,1+(j-1)/stepi)=NaN;
   else
    X(1+(i-1)/stepi,1+(j-1)/stepi)=mean(submat(:));
   end
  end; end
  x=X;

  for i=1:stepi:size(y,1); for j=1:stepi:size(y,2)
   submat=y(i:min(i+stepi-1,size(y,1)),j:min(j+stepi-1,size(y,2)));
   if any(isnan(submat(:)))
    Y(1+(i-1)/stepi,1+(j-1)/stepi)=NaN;
   else
    Y(1+(i-1)/stepi,1+(j-1)/stepi)=mean(submat(:));
   end
  end; end
  y=Y;

  for i=1:stepi:size(u,1); for j=1:stepi:size(u,2)
   submat=u(i:min(i+stepi-1,size(u,1)),j:min(j+stepi-1,size(u,2)));
   if sum(isnan(submat(:)))>=nsamples/2
    U(1+(i-1)/stepi,1+(j-1)/stepi)=NaN;
   else
    U(1+(i-1)/stepi,1+(j-1)/stepi)=median(submat(~isnan(submat(:))));
   end
  end; end
  u=U;

  for i=1:stepi:size(v,1); for j=1:stepi:size(v,2)
   submat=v(i:min(i+stepi-1,size(v,1)),j:min(j+stepi-1,size(v,2)));
   if sum(isnan(submat(:)))>=nsamples/2
    V(1+(i-1)/stepi,1+(j-1)/stepi)=NaN;
   else
    V(1+(i-1)/stepi,1+(j-1)/stepi)=median(submat(~isnan(submat(:))));
   end
  end; end
  v=V;

 % or thin vectors by sub-sampling or taking mean
 else

  if method==2
    disp(['Thinning vectors with means of ',num2str(nsamples),' samples per vector'])
  elseif method==1
    disp(['Thinning vectors by sub-sampling one from ',num2str(nsamples),' vectors'])
  else
    error('method has an incorrect value')
  end

  for k=1:max(1,times)

   X=ones(floor(size(x,1)/2), floor(size(x,2)/2));
   for i=1:2:size(x,1)-1 ; for j=1:2:size(x,2)-1 ;
	if method==1
	 X((i+1)/2,(j+1)/2)=x(i,j);
	else
	 X((i+1)/2,(j+1)/2)=mean2(x(i:i+1,j:j+1));
	end
   end; end;
   x=X;
   Y=ones(floor(size(y,1)/2), floor(size(y,2)/2));
   for i=1:2:size(y,1)-1 ; for j=1:2:size(y,2)-1 ;
	if method==1
 	 Y((i+1)/2,(j+1)/2)=y(i,j);
	else
 	 Y((i+1)/2,(j+1)/2)=mean2(y(i:i+1,j:j+1));
	end
   end; end;
   y=Y;
   U=ones(floor(size(u,1)/2), floor(size(u,2)/2));
   for i=1:2:size(u,1)-1 ; for j=1:2:size(u,2)-1 ;
	if method==1
	 U((i+1)/2,(j+1)/2)=u(i,j);
	else
	 U((i+1)/2,(j+1)/2)=mean2(u(i:i+1,j:j+1));
	end
   end; end;
   u=U;
   V=ones(floor(size(v,1)/2), floor(size(v,2)/2));
   for i=1:2:size(v,1)-1 ; for j=1:2:size(v,2)-1 ;
	if method==1
	 V((i+1)/2,(j+1)/2)=v(i,j);
	else
	 V((i+1)/2,(j+1)/2)=mean2(v(i:i+1,j:j+1));
	end
   end; end;
   v=V;

  end

 end

end

if ismap(h)

 [lat,lon] = minvtran(x,y);

 x=lat;
 y=lon;

end

if debug; disp(['...Leaving ',mfilename]); end

