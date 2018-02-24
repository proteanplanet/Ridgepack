function ridgepack_mask(a,b,zi,col,ref,lw)

% ridgepack_mask - Adds a black border around NaN values on a plot
%
% function ridgepack_mask(a,b,zi,col,ref,lw)
% 
% This function draws a border on a Cartesian plot around regions of 
% NaN or a specific reference value.  It is specifically intended
% to highlight the land-sea boundaries or specific values in a field,
% such as a zero dividing line, in numerical model output from 
% climate, weather and ocean models.
%
% Note that sometimes the mask appears to not align with 
% a mask that is filled, using ridgepack_image, but this is typically
% a matlab graphics problem, not a misallignment.
%
% INPUT:
%
% a   - x-direction coordinates
% b   - y-direction coordinates
% zi  - values of axb array containing NaN or reference values as a mask
% col - optional input for color of dividing line (e.g. 'k' for black)
% ref - reference value if it is other than NaN (default is NaN)
% lw  - line width (default is 0.4)
%
% Ridgepack Version 1.0
% Andrew Roberts, Naval Postgraduate School, March 2018 (afrobert@nps.edu)


global debug;
if debug; disp(['Entering ',mfilename,'...']); end

if nargin<3; error('Not enough arguments'); end
if nargin<4; col='k'; end
if nargin<5; ref=NaN ; end
if nargin<6; lw=0.4 ; end

if not(ndims(zi)==2);error('zi is not a 2D array');end

if length(a)~=size(zi,1) || length(b)~=size(zi,2); 
 disp(['length of a is ',num2str(length(a))])
 disp(['length of b is ',num2str(length(b))])
 disp(['size of zi is ',num2str(size(zi))])
 error('size of axes and data incorrect');
end 

if debug; disp([mfilename,'1']); end

% extend the arrays to account for boundaries
ii=size(zi,1);
jj=size(zi,2);

z(1:ii+2,1:jj+2)=NaN;
z(2:ii+1,2:jj+1)=zi(:,:);

x(2:ii+1)=a(:); x(1)=a(1); x(ii+2)=a(ii);
y(2:jj+1)=b(:); y(1)=b(1); y(jj+2)=b(jj);

% assemble coordinates of land edge
xc=(x(1:ii+1)+x(2:ii+2))/2;
yc=(y(1:jj+1)+y(2:jj+2))/2;

% create array tracing the mask
if isnan(ref)
 xxor(2:ii+1,2:jj+1)=xor(isnan(z(2:ii+1,2:jj+1)),isnan(z(3:ii+2,2:jj+1)));
 yxor(2:ii+1,2:jj+1)=xor(isnan(z(2:ii+1,2:jj+1)),isnan(z(2:ii+1,3:jj+2)));
else
 xxor(2:ii+1,2:jj+1)=xor(z(2:ii+1,2:jj+1)<ref,z(3:ii+2,2:jj+1)<ref);
 yxor(2:ii+1,2:jj+1)=xor(z(2:ii+1,2:jj+1)<ref,z(2:ii+1,3:jj+2)<ref);
end

% create mask
cxor=ones(size(xxor));
cyor=ones(size(yxor));
cxor(not(xxor))=NaN;
cyor(not(yxor))=NaN;

if debug; disp([mfilename,'2']); end

% generate lines
xx=zeros((ii+1)*(jj+1)*5,1);
yy=zeros((ii+1)*(jj+1)*5,1);
k=1;
for i=2:ii+1;
for j=2:jj+1;

  xx(k)=cyor(i,j)*xc(i-1); yy(k)=cyor(i,j)*yc(j);
  xx(k+1)=cyor(i,j)*xc(i); yy(k+1)=cyor(i,j)*yc(j);

  xx(k+2)=cxor(i,j)*xc(i); yy(k+2)=cxor(i,j)*yc(j);
  xx(k+3)=cxor(i,j)*xc(i); yy(k+3)=cxor(i,j)*yc(j-1);

  xx(k+4)=NaN; yy(k+4)=NaN;

  k=k+5;

end
end

line(xx,yy,'Color',col,'linewidth',lw);

drawnow

if debug; disp(['...Leaving ',mfilename]); end

