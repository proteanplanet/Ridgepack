function [cont]=ridgepack_contlev(z,v)

% ridgepack_contlev - Generates contour levels automatically
%
% function [cont]=ridgepack_contlev(z,v)
%
% INPUT:
% 
% z    - 2D matrix of values to be contoured
% v    - vector of contours to be calculated (optional)
%
%
% OUTPUT:
%
% cont - matrix of contour levels
%
% Note that this function uses the matlab contourc function
% to generate the contour levels, and this can, on occasions
% produce some unexpected results near zero. For this reason
% contours are rounded to the nearest six smalled orders of 
% magnitude from the mean contour value.
%
% Ridgepack Version 1.0
% Andrew Roberts, Naval Postgraduate School, March 2018 (afrobert@nps.edu)
% Blair Greenan's code segment added from Matlab Central

global debug;
if debug; disp(['Entering ',mfilename,'...']); end

% generate contour levels
if nargin<2
 cout=contourc(double(z));
else
 if length(v)<2; error('Contour vector must have at least 2 elements'); end
 cout=contourc(double(z),v);
end

if isempty(cout)
	if isempty(z)
		error('z is empty including when masked out')
	else
		error('cout is empty - z is a constant field including when masked out')
	end
end

% get contour values automatically if cont not suppplied 
% (code copied from Blair Greenan's colorbarf.m at matlab central)
i = 1;
while ~isempty(cout)
   C1(i) = cout(1,1); % contour level
   N1 = cout(2,1); % number of points in contour
   cout(:,1:N1+1) = []; % shrink the matrix
   if (N1 > 1) 
     i= i + 1;
   end
end
C1=C1(:);
cont(1)=C1(1);
j=1;
for i=2:length(C1)
 if C1(i)~=C1(i-1)
  j=j+1;
  cont(j)=C1(i);
 end
end 

% split contour levels to double number of bands
oldcont=cont;
j=1;
for i=2:length(oldcont)
 j=j+1;
 cont(j)=(oldcont(i)+oldcont(i-1))/2;
 j=j+1;
 cont(j)=oldcont(i);
end

% round contours to the nearest 6 smaller orders of magnitude 
% to make allowance for problems close to zero for contourc.
cont=round((10^6)*cont)/(10^6);

% make sure there are at least two contour levels selected
if length(cont)<2; error('Less than two contour levels calculated'); end

if debug; disp(['...Leaving ',mfilename]); end
