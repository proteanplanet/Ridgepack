function [cmap]=ridgepack_colormap(cont,ref,colors,logscale)

% function [cmap]=ridgepack_colormap(cont,ref,colors,logscale)
%
% This function is part of Ridgepack Version 1.0.
% It creates the colormap and draws a colorbar to the right
% of or underneath the plot. Units are added, and a factor of ten
% is also added if the number is less than 1 or greater than or
% equal to 1000.  The colormap and colorbar have several features
% that improve the style and readability of the plots, and ability
% to distinguish between values, which can be difficult using 
% standard matlab color scales. Note that it is assumed that the 
% data plotted is already indexed data which can be obtained
% by running the nccolorindex function on data before 
% using them in color utilities such as pcolor.
%
% INPUTS:
%
% cont   - color levels entered as a vector [C1,C2,...,CX].
%
% ref    - reference point about which the color is centered, and, when
%          the minimum or maximum is equal to ref, all points below or 
%          above this value, respectively, are not filled with a color.
%          {zero is the default value}
%
% colors - This sets the color scheme required:
%	   'bluered'   - highest contrast blue fading through red at ref value
%          'greenred'  - green fading through red at ref value
%          'bluegreen' - blue fading through green at ref value
%          'jet'       - standard matlab 'jet' color scheme.  Ref has no effect.
%          'cool'      - standard matlab 'cool' color scheme.  Ref has no effect.
%          'gray'      - straight gray color scale. Ref has no effect.
%          'parula'    - standard matlab 'parula' color scheme. Ref has no effect.
%          'summer'    - standard matlab 'summer' color scheme. Ref has no effect.
%          {bluered is the default}
%
% logscale - set to true if log scale is being used (optional).
%
% OUTPUT:
%
% cmap   - colormap output
%
% The only two inputs that are compulsory are cont and units, the
% others can be eliminated if need be.
% 
% Andrew Roberts, Naval Postgraduate School, March 2018 (afrobert@nps.edu)

ht=get(gcf,'CurrentAxes');

% change renderer
set(gcf,'renderer','zbuffer')

% set defaults
if nargin<1; error('cont not specified'); end
if nargin<2; ref=0.0; end
if nargin<3; colors='bluered' ; end
if nargin<4; logscale=false ; end

% sort cont in ascending numerical order, taking into account a single value
if isempty(cont)
  error('Cont is an empty vector');
elseif length(cont)<2
  disp('Contour levels should be at least two elements long'); 
  c1=cont(1); 
  cont=zeros(1, 3);
  cont(1)=c1-eps;
  cont(2)=c1;
  cont(3)=c1+eps;
end
cont=sort(cont);

% set tick marks
Ytick=zeros(1, length(cont));

% set new color schemes
cmap=1:1:length(cont)-1;
zp=ceil(median(1:1:length(cmap)));
yp=zp;

for i=2:length(cmap);
 if cont(i)>ref & cont(i-1)<ref
   if logscale
    zp=i-1;  
   else
    zp=i;    
   end
   yp=i-1;
 elseif cont(i)>=ref & cont(i-1)<ref
   zp=i;
   yp=i;
 elseif cont(i)<=ref & cont(i-1)>ref
   error('contour color scale is reversed - not allowed');
 end
end


% build colormap 
if strcmp(colors,'bluered');

 % set initial dummy colormap
 colmap=colormap(jet(length(cmap))); 

 if size(colmap,1)~=length(cmap) ; 
  error('cmap programming error'); 
 end

 huem=rgb2hsv(colmap);
 vb=(length(cmap)-1)/(max(yp,length(cmap)-zp)); 

 for i=1:length(cmap);
  if i>1 & i>=zp;
   coef=(i-zp)/(length(cmap)-1);
   huem(i,1)=0.20-0.30*vb*coef;
   huem(i,2)=0.15+0.85*vb*coef;
   huem(i,3)=1.0-0.5*vb*coef;
  elseif i<yp
   coef=(yp-i-1)/(length(cmap)-1);
   huem(i,1)=0.45+0.30*vb*coef;
   huem(i,2)=0.15+0.85*vb*coef;
   huem(i,3)=1.0-0.5*vb*coef;
  end
  huem(i,1)=max(min(huem(i,1),1),0);
  huem(i,2)=max(min(huem(i,2),1),0);
  huem(i,3)=max(min(huem(i,3),1),0);
 end

 cmap=hsv2rgb(huem);
 colormap(cmap);

elseif strcmp(colors,'greenred');

 % set initial dummy colormap
 colmap=colormap(jet(length(cmap))); 

 if size(colmap,1)~=length(cmap) ; 
  error('cmap programming error'); 
 end

 huem=rgb2hsv(colmap);
 vb=(length(cmap)-1)/(max(yp,length(cmap)-zp)); 

 for i=1:length(cmap);
  if i>1 & i>=zp;
   coef=(i-zp)/(length(cmap)-1);
   huem(i,1)=0.20-0.30*vb*coef;
   huem(i,2)=0.15+0.80*vb*coef;
   huem(i,3)=1.0-0.5*vb*coef;
  elseif i<yp
   coef=(yp-i-1)/(length(cmap)-1);
   huem(i,1)=0.30+0.30*vb*coef;
   huem(i,2)=0.15+0.80*vb*coef;
   huem(i,3)=1.0-0.5*vb*coef;
  end
  huem(i,1)=max(min(huem(i,1),1),0);
  huem(i,2)=max(min(huem(i,2),1),0);
  huem(i,3)=max(min(huem(i,3),1),0);
 end

 cmap=hsv2rgb(huem);
 colormap(cmap);

elseif strcmp(colors,'bluegreen');

 % set initial dummy colormap
 colmap=colormap(jet(length(cmap))); 

 if size(colmap,1)~=length(cmap) ; 
  error('cmap programming error'); 
 end

 huem=rgb2hsv(colmap);
 vb=(length(cmap)-1)/(max(yp,length(cmap)-zp)); 

 for i=1:length(cmap);
  if i>1 & i>=zp;
   coef=(i-zp)/(length(cmap)-1);
   huem(i,1)=0.45-0.20*vb*coef;
   huem(i,2)=0.15+0.80*vb*coef;
   huem(i,3)=1.0-0.5*vb*coef;
  elseif i<yp
   coef=(yp-i-1)/(length(cmap)-1);
   huem(i,1)=0.55+0.20*vb*coef;
   huem(i,2)=0.15+0.80*vb*coef;
   huem(i,3)=1.0-0.5*vb*coef;
  end
  huem(i,1)=max(min(huem(i,1),1),0);
  huem(i,2)=max(min(huem(i,2),1),0);
  huem(i,3)=max(min(huem(i,3),1),0);
 end

 cmap=hsv2rgb(huem);
 colormap(cmap);

elseif strcmp(colors,'jet'); % standard matlab jet color scheme

 cmap=colormap(jet(length(cmap))); 

elseif strcmp(colors,'cool'); % standard matlab jet color scheme

 cmap=colormap(cool(length(cmap))); 

elseif strcmp(colors,'parula'); % standard matlab jet color scheme

 cmap=colormap(parula(length(cmap))); 

elseif strcmp(colors,'summer'); % standard matlab jet color scheme

 cmap=colormap(summer(length(cmap))); 

else %  grey colormap

 colormap(contrast(length(cmap)))
 brighten(0.5)

end


% Change color to grey over the zero or reference
% contour if it is not exactly specified, otherwise
% just draw a white line along zero or reference
if zp~=yp && any(strcmpi(colors,{'bluered','greenred','bluegreen','jet','cool','parula','summer'})); 
 for i=1:length(Ytick);
  if zp==i
   cmap(i-1,:)=0.90;
  end
 end
 colormap(cmap);
end


% get final colormap
cmap=colormap;

