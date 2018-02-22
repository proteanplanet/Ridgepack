function ht=ridgepack_text(x,y,string,fontsize,color,halign,valign,backgroundcolor,interpreter)

% ridgepack_text - Text plotting function to annotate figures along lines
%
% function ht=ridgepack_text(x,y,string,fontsize,color,halign,valign,backgroundcolor,interpreter)
%
% This function plots curved text strings along given lines on plots specified 
% with the Cartesian coordinates x and y. This can be used to mark
% contours, or add annotations along lines on graphs.  The difference
% between this and the standard MATLAB text function is that the characters
% align with curves provide, rather than plotting text in a straight line.  
% x and y must be vectors, rather than a single set of coordinates. If not,
% then the regular MATLAB text file should be used.  This function works with
% both linear and logarithmic axes. When plotting equations or symbols
% in using the default Latex interpreter, text buttressed by math mode indicators
% (i.e. '$$' or '$') will be plotted at a constant angle.
%
% INPUT:
%
% x               - x-coordinate (Cartesian) 
% y               - y-coordinate (Cartesian)
% string          - Text string for annotation
% fontsize        - Fontsize in points (optional; default is axes fontsize)
% color           - Color of text specified either as a RGB triplet
%                   e.g. [0 0 0](optional; default is black)
% halign          - Horizontal text alignment (left, center, right)
%                   (optional; default is center)
% valign          - Vertical text alignment (top, middle, bottom, cap, baseline)
%                   (optional; default is middle)
% backgroundcolor - Color of text box specified either as a RGB triplet
%                   e.g. [1 1 1], or 'none' (optional; default is [1 1 1])
% interpreter     - Text interpreter (optional; default is 'Latex')
%
%
% OUTPUT:
%
% ht - Cell array of handles for each character being plotted.
%
% IMPORTANT:
% This function must be called after the axes aspect ratios and limits
% have been established.  It does not have a callback function for realigning
% text if changes are made subsequent to plotting text.
%
% Notes:
% In future, this function could easily be expanded to work on maps for
% x and y representing latitude and longitude.
% 
% Ridgepack Version 1.0
% Andrew Roberts, Naval Postgraduate School, March 2018 (afrobert@nps.edu)

% check inputs
if nargin<3
 error('insufficient arguments')
end
ha=get(gcf,'CurrentAxes');
if nargin<4
 if isempty(ha)
  axes;
 end
 fontsize=get(gca,'FontSize');
end
if nargin<5
 color=0*[1 1 1];
end
if nargin<6
 halign='center';
end
if nargin<7
 valign='middle';
end
if nargin<8
 backgroundcolor=[1 1 1];
end
if nargin<9
 interpreter='Latex';
end

% convert coordinates if using a map
if ismap(ha)

 disp('Adding text to a map')
 mstruct=gcm;
 if strmatch(mstruct.mapprojection,'globe')
   error('This function does not work with the globe projection')
 end
 [x,y]=mfwdtran(mstruct,x,y,ha,'surface'); % untrimmed values

end

% build high resolution x and y interpolated functions
if length(x)>1 & length(x)==length(y)
 idx=0;
 for i=1:length(x)-1
  for j=1:1000;
   idx=idx+1;
   xx(idx)=x(i)+(j-1)*(x(i+1)-x(i))./1000;
   yy(idx)=y(i)+(j-1)*(y(i+1)-y(i))./1000;
  end
 end
 x=xx;
 y=yy;
% s=sqrt((x(2:end)-x(1:end-1))^2+(y(2:end)-y(1:end-1))^2);
 
elseif length(x)~=length(y)
 error('length of x and y is different')
elseif length(x)<2
 error('length of x and y is 1, use normal MATLAB text function')
end

% get plot box aspect ratio
pos=get(gca,'Position');
pbar=get(gca,'PlotBoxAspectRatio');
xlims=xlim;
ylims=ylim;
if strcmp(get(gca,'XScale'),'log')
 dar(1)=diff(log10(xlims));
else
 dar(1)=diff(xlims);
end
if strcmp(get(gca,'YScale'),'log')
 dar(2)=diff(log10(ylims));
else
 dar(2)=diff(ylims);
end
if strcmp(get(gca,'XDir'),'reverse')
 dar(1)=-dar(1);
end
if strcmp(get(gca,'YDir'),'reverse')
 dar(2)=-dar(2);
end
coeff=(pbar(2)/pbar(1))*(dar(1)/dar(2));

% smooth the contour data
if strcmp(get(gca,'XScale'),'log')
 xrough=log10(x);
else
 xrough=x;
end
if strcmp(get(gca,'YScale'),'log')
 yrough=log10(y);
else
 yrough=y;
end
xsmooth=[xrough(1):(xrough(end)-xrough(1))/1000:xrough(end)];
ysmooth=pchip(xrough,yrough,xsmooth);
if strcmp(get(gca,'XScale'),'log')
 x=(10.^xsmooth);
else
 x=xsmooth;
end
if strcmp(get(gca,'YScale'),'log')
 y=(10.^ysmooth);
else
 y=ysmooth;
end

% estimate starting location
if strcmp(halign,'left')
 idx=3;
elseif strcmp(halign,'center')
 if strcmp(get(gca,'XScale'),'log') & strcmp(get(gca,'YScale'),'log')
  anglet=atan2d((log10(y(floor(end/2)+400))-log10(y(floor(end/2)-400)))*coeff,...
                 log10(x(floor(end/2)+400))-log10(x(floor(end/2)-400)));
 elseif strcmp(get(gca,'XScale'),'log') 
  anglet=atan2d((y(floor(end/2)+400)-y(floor(end/2)-400))*coeff,...
                 log10(x(floor(end/2)+400))-log10(x(floor(end/2)-400)));
 elseif strcmp(get(gca,'YScale'),'log') 
  anglet=atan2d((log10(y(floor(end/2)+400))-log10(y(floor(end/2)-400)))*coeff,...
                 x(floor(end/2)+400)-x(floor(end/2)-400));
 else
  anglet=atan2d((y(floor(end/2+400))-y(floor(end/2-400)))*coeff,...
                 x(floor(end/2+400))-x(floor(end/2-400)));
 end
 h=text(x(floor(end/2)),y(floor(end/2)),string,...
                       'Rotation',anglet,...
                       'Color',color,...
                       'BackgroundColor',backgroundcolor,...
                       'Margin',1,...
                       'Units','data',...
                       'FontSize',fontsize,...
                       'HorizontalAlignment','center',...
                       'VerticalAlignment',valign,...
                       'Interpreter',interpreter);
 ext=get(h,'Extent');
 delete(h);
 if anglet>45 & anglet<135
  idx=find(abs(y-ext(2))==min(abs(y-ext(2))));
 else
  idx=find(abs(x-ext(1))==min(abs(x-ext(1))));
 end
elseif strcmp(halign,'right')
 if strcmp(get(gca,'XScale'),'log') & strcmp(get(gca,'YScale'),'log')
  anglet=atan2d((log10(y(end))-log10(y(end-400)))*coeff,...
                 log10(x(end))-log10(x(end-400)));
 elseif strcmp(get(gca,'XScale'),'log')
  anglet=atan2d((y(end)-y(end-400))*coeff,log10(x(end))-log10(x(end-400)));
 elseif strcmp(get(gca,'YScale'),'log')
  anglet=atan2d((log10(y(end))-log10(y(end-400)))*coeff,x(end)-x(end-400));
 else
  anglet=atan2d((y(end)-y(end-400))*coeff,x(end)-x(end-400));
 end

 h=text(x(end),y(end),string,...
                       'Rotation',anglet,...                       
                       'Color',color,...
                       'Margin',0.00001,...
                       'Units','data',...
                       'FontSize',fontsize,...
                       'HorizontalAlignment','right',...
                       'VerticalAlignment',valign,...
                       'Interpreter',interpreter);
 ext=get(h,'Extent');
 delete(h);
 if anglet>45 & anglet<135
  idx=find(abs(y-ext(2))==min(abs(y-ext(2))));
 else
  idx=find(abs(x-ext(1))==min(abs(x-ext(1))));
 end
else
 error('can only choose left, right or center location for text')
end

% determine blocks of text within Latex equation descriptors 
index=1;
block=0;
while index<=length(string)

 k=(index-1)+strfind(string(index:end),'$$');
 if ~isempty(k) & length(k)>1 & k(1)==index
  block=block+1;
  stringblock{block}=string(k(1):k(2)+1);
  colorblock{block}=color;
  index=k(2)+2;
 elseif ~isempty(k) & length(k)==1 
  error('Found a set of unended Latex $$ descriptor in string')
 else
  k=(index-1)+strfind(string(index:end),'$');
  if ~isempty(k) & length(k)>1 & k(1)==index
   block=block+1;
   stringblock{block}=string(k(1):k(2));
   colorblock{block}=color;
   index=k(2)+1;
  elseif ~isempty(k) & length(k)==1
   error('Found an unended Latex $ descriptor in string')
  else
   block=block+1;
   if strcmp(string(index),' ') & index>1 & index<length(string) & ...
     ~(strcmp(string(index-1),'$') & strcmp(string(index+1),'$'))
    stringblock{block}='.';
    colorblock{block}=backgroundcolor;
    index=index+1;
   else
    stringblock{block}=string(index);
    colorblock{block}=color;
    index=index+1;
   end
  end
 end

end

% plot text one character at a time from left start
deleteall=false;
anglet=0;
idxmin=idx;
for i=1:block

  % find the nearest point on the high resolution grid
  if i>1
   if anglet>45 & anglet<135
    idx=find(abs(y-(ext(2)+ext(3)*sind(anglet)))==...
             min(abs(y-(ext(2)+ext(3)*sind(anglet)))));
   else
    idx=find(abs(x-(ext(1)+ext(3)*cosd(anglet)))==...
             min(abs(x-(ext(1)+ext(3)*cosd(anglet)))));
   end
  end

  % find text extent details
  h=text(x(idx),y(idx),char(stringblock{i}),...
                       'Rotation',0,...
                       'Color',color,...
                       'Margin',1,...
                       'Units','data',...
                       'FontSize',fontsize,...
                       'HorizontalAlignment','left',...
                       'VerticalAlignment',valign,...
                       'Interpreter',interpreter);
  ext=get(h,'Extent');
  delete(h);

  % find index of text on x-y grid for angled text
  idxangle=idx+20*length(char(stringblock{i}));
  idxangle=min(idxangle,length(x));

  if idxangle==idx
   disp(['WARNING: Unable to fit text to the given coordinate space: ',string])
  elseif strcmp(get(gca,'XScale'),'log') & strcmp(get(gca,'YScale'),'log')
   anglet=atan2d((log10(y(idxangle))-log10(y(idx)))*coeff,...
                  log10(x(idxangle))-log10(x(idx)));
  elseif strcmp(get(gca,'XScale'),'log')
   anglet=atan2d((y(idxangle)-y(idx))*coeff,log10(x(idxangle))-log10(x(idx)));
  elseif strcmp(get(gca,'YScale'),'log')
   anglet=atan2d((log10(y(idxangle))-log10(y(idx)))*coeff,x(idxangle)-x(idx));
  else
   anglet=atan2d((y(idxangle)-y(idx))*coeff,x(idxangle)-x(idx));
  end

  % plot white underlay over lines to avoid breaks between text
  % (this is done at a smaller font size to avoid end blowouts)
  if i==1 & ~ischar(backgroundcolor)
   gt=text(x(idx),y(idx),string,...
                       'Rotation',anglet,...
                       'Color',backgroundcolor,...
                       'BackgroundColor',backgroundcolor,...
                       'Margin',0.0001,...
                       'Units','data',...
                       'FontSize',0.9*fontsize,...
                       'EdgeColor','none',...
                       'HorizontalAlignment','left',...
                       'VerticalAlignment',valign,...
                       'Interpreter',interpreter);
  else
   gt=[];
  end
  
  % plot text
  ht{i}=text(x(idx),y(idx),char(stringblock{i}),...
                       'Rotation',anglet,...
                       'Color',colorblock{i},...
                       'BackgroundColor',backgroundcolor,...
                       'Margin',0.0001,...
                       'Units','data',...
                       'FontSize',fontsize,...
                       'EdgeColor','none',...
                       'HorizontalAlignment','left',...
                       'VerticalAlignment',valign,...
                       'Interpreter',interpreter);

  % make sure all text is within the plotting frame
  extcheck=get(ht{i},'Extent');
  if extcheck(1)<=xlims(1) | extcheck(2)<=ylims(1) | ...
     extcheck(1)+extcheck(3)>=xlims(2) | extcheck(2)+extcheck(4)>=ylims(2)
   deleteall=true;
   break
  end

end

% delete strings not shown
if deleteall
 disp(['REMOVED: String ',string,' falls at least partly outside of plot box'])
 if ishandle(gt)
  delete(gt)
 end
 for j=1:i
  delete(ht{j})
 end
% or add blank string behind
elseif ishandle(gt)
 ht{i+1}=gt;
end

if nargout==0
 clear ht
end

