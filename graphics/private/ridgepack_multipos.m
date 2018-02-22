function [ha]=ridgepack_multipos(hf,ha,nrows,ncols,row,col,notation,fontsize,fontcol,fontbox)

% ridgepack_multipos - Main positioning engine for constructing multiple axes figures
%
% function [ha]=ridgepack_multipos(hf,ha,nrows,ncols,row,col,notation,fontsize,fontcol,fontbox)
%
% This function positions axes according to the row and column in a multiple
% axes figure. It also adds two extra, invisible axes, to the figure, used
% for positioning a global colorbar, and a title header. This function
% is called by ridgepack_multiplot and ridgepack_multialign, and is not available directly
% for the user to call. 
%
% Inputs:
% hf       - Figure handle
% ha       - Axes handle (may be left empty if no axes yet created)
% nrows    - Number of rows of axes in the figure
% ncols    - Number of columns of axes in the figure
% row      - Row of current axes being positioned
% col      - Column of current axes being positioned
% notation - Text with notation within each subplot. This is a char variables,
%            but when supplied as a number, indicates to plot the notation
%            that is stored away until calling ridgepack_multialign (parent function)
% fontsize - Manually set font size for the notation (typically 8-14; optional)
% fontcol  - Manually set font color for the notation (RBG vector; optional)
% fontbox  - If set to 1, provides a box around the notation, if set to 2,
%            puts notation in lower right, if set to 3, puts it in the lower
%            left, and if set to 4, puts it in the upper right. If set to 5
%            it places the number 1/4 down the left side.
%
%
% OUTPUT:
%
% ha    - The axis handle of the row,col axes used or created during execution
%
% If ncols==1 (a single column figure), then the default span of the axes
% is to stretch across the entire figure, as is useful for comparing graphs
% stacked upon each other.
%
% This function, aside from positioning axes, assigns a data structure
% that describes the multiplot configuration to app data for the given
% figure window.
%
% Ridgepack Version 1.0
% Andrew Roberts, Naval Postgraduate School, March 2018 (afrobert@nps.edu)

global debug; if isempty(debug); debug=false; end
if debug; disp(['Entering ',mfilename,'...']); end

% default
plotnotation=false;

% briefly check inputs (already checked in ridgepack_multiplot)
if nargin<6
 error('incorrect number of inputs')
elseif nargin>=7 & isnumeric(notation)
 plotnotation=true;
elseif nargin>=7 & ~ischar(notation)
 error('notation must be a string')
end

% check font size and font color
if nargin>7 & ~isnumeric(fontsize)
 error('Font size for notation must be a number')
elseif nargin>8 & ~isnumeric(fontcol)
 error('Font color for notation must be a number')
elseif nargin>8 & length(fontcol)~=3
 error('Font color must be a three element RGB vector')
elseif nargin<=7
 fontsize=12;
 fontcol=[0 0 0];
elseif nargin<=8
 fontcol=[0 0 0];
end

% more error checking
if length(nrows)>1 | length(ncols)>1
 error('nrows and ncols must be a single number, not a matrix or vector')
end

% add fontbox
if nargin>=10 & ~isnumeric(fontbox)
 error('fontbox is not a number')
elseif nargin<10
 fontbox=0;
end

% get set minumum gap between plots and edge in normalized coordinates
mingap=0.005/min(2,sqrt(max(1,max(nrows,ncols)-1)));
minedge=0.01/min(2,sqrt(max(1,max(nrows,ncols)-1)));

% set units to use to create multi plot figure (normalized)
units='normalized';

% global colorbar settings

% initial setting of colorbar span relative to figure
cbspan=0.90; 

% minumum separation from axes for colorbar
% (must be same as in private/ridgepack_cbpos.m)
minsep=0.01/min(2,sqrt(max(1,max(nrows,ncols)-1))); 

if ~ishandle(hf) 
 error('hf is not a handle')
elseif ~isempty(ha) & ~ishandle(ha)
 error('ha is not a handle')
end

% get extent of figure in pixels
hfpos=getpixelposition(hf);
widthfigpixels=hfpos(3); 
heightfigpixels=hfpos(4); 

% get width (rx) and height (ty) of plot
if strcmpi(units,'normalized')
 rx=1; 
 ty=1;
elseif strcmpi(units,'pixels')
 rx=widthfigpixels;
 ty=heightfigpixels;
 mingap=mingap*max(rx,ty);
 minedge=minedge*max(rx,ty);
else
 error('Units not recognized')
end

% set maximum of rows and columns
maxrc=max(nrows,ncols);

% get configuration
r=getappdata(hf,'MultiplotConfiguration');

% set or reset based configuration if required
if isempty(r) 

 % assign default extents
 if ncols==1
  extentx=rx-2*minedge-2*mingap;
 else
  extentx=(rx-(2*minedge)-(mingap*(maxrc-1)))/maxrc;
 end
 extenty=(ty-(2*minedge)-(mingap*(maxrc-1)))/maxrc;
 r.extent.x=extentx*ones(ncols,1);
 r.extent.y=extenty*ones(nrows,1);

 % assign dead space fillers
 r.dead.x=zeros(nrows,ncols);
 r.dead.y=zeros(nrows,ncols);

 % assign defaults to the multiplot configuration structure
 r.tightinset=mingap*ones(nrows,ncols,4);
 r.tightinset(:,1,1)=minedge;
 r.tightinset(:,end,3)=minedge;
 r.tightinset(end,:,2)=minedge;
 r.tightinset(1,:,4)=minedge;

 % assign default gaps
 r.gap.x=mingap*ones(ncols+1,1);
 r.gap.y=mingap*ones(nrows+1,1);

 % assign dummy handles and 
 r.handle=NaN*zeros(nrows,ncols); 

 % assign colorbar flags
 r.colorbar=zeros(nrows,ncols);

 % assign notation text
 r.notation=cell(nrows,ncols);
 r.notationhandle=NaN*zeros(nrows,ncols);

 % create new multiplot axes for the whole plot for the colorbar
 r.multihandle=axes('Position',[minedge/rx minedge/ty 1-2*minedge/rx 1-2*minedge/ty],...
          'Visible','off','XTickLabel',[],'YTickLabel',[],'Box','on','Layer','top');

 % create an axis for the main title and legend placement
 r.titlehandle=axes('Position',[minedge/rx minedge/ty 1-2*minedge/rx 1-2*minedge/ty],...
          'Visible','on','XTick',[],'YTick',[],'XTickLabel',[],'YTickLabel',[],...
          'Color','none','Box','on','XColor',[1 1 1],'YColor',[1 1 1]);

 % turn off the title axes unless debugging
 if ~debug ; axis(r.titlehandle,'off'); end

 % assign MATLAB legend and colorbar handle field
 r.legcolhandle=NaN*zeros(nrows,ncols);

 % write app data
 setappdata(hf,'MultiplotConfiguration',r)

elseif ~isstruct(r)
 error('structure is not assigned to MultiplotConfiguration')
elseif ~isfield(r,'gap')
 error('gap is not a field in the MultiplotConfiguration')
elseif ~isfield(r.gap,'x') | ~isfield(r.gap,'y')
 error('gap x or y are not fields in the MultiplotConfiguration')
elseif ~isfield(r,'extent')
 error('extent is not a field in the MultiplotConfiguration')
elseif ~isfield(r.extent,'x') | ~isfield(r.extent,'y')
 error('extent x or y are not fields in the MultiplotConfiguration')
elseif length(r.extent.x)~=ncols | length(r.extent.y)~=nrows
 error('MultiplotConfiguration ncols or nrows inconsistent')
elseif ~isfield(r,'tightinset')
 error('tightinset is not a field in the MultiplotConfiguration')
elseif size(r.tightinset)~=[nrows ncols 4]
 error('tightinset is not a field in the MultiplotConfiguration')
elseif ~isfield(r,'handle')
 error('handle is not a field in the MultiplotConfiguration')
elseif ~isfield(r,'colorbar')
 error('colorbar is not a field in the MultiplotConfiguration')
elseif ~isfield(r,'multihandle')
 error('handle of multiplot global axes not found')
elseif ~isfield(r,'notation')
 error('notation field not found')
elseif ~isfield(r,'notationhandle')
 error('notation handle field not found')
elseif ~ishandle(r.multihandle)
 error('multiplot colorbar axes handle error ')
elseif ~isfield(r,'titlehandle')
 error('handle of multiplot title axes not found')
elseif ~ishandle(r.titlehandle)
 error('multiplot title axes handle error ')
elseif ~isfield(r,'legcolhandle')
 error('no MATLAB legend and colorbar handle field found')
end

% remove existing axes if creating a new one
if isempty(ha) & ~isnan(r.handle(row,col)); delete(ha); end

% set units to those required
if ishandle(ha); haunits=get(ha,'Units'); set(ha,'Units',units); end

% turn on box around ploit
if ishandle(ha); set(ha,'box','on'); end

% Get tight inset information for the current axes
if ishandle(ha); r=ridgepack_multitight(ha,r,units,mingap,row,col); end

% Get tight inset information from previously set axes
for ro=1:nrows; for co=1:ncols
 if ishandle(r.handle(ro,co)); r=ridgepack_multitight(r.handle(ro,co),r,units,mingap,ro,co); end
end; end

% get the space required for the main title
tunits=get(r.titlehandle,'Units');
set(r.titlehandle,'Units',units)
titleinset=get(r.titlehandle,'TightInset');

% determine if the global axes has a colorbar and its required edge gap
gunits=get(r.multihandle,'Units');
set(r.multihandle,'Units',units);
gposition=get(r.multihandle,'Position');
gtightinset=get(r.multihandle,'TightInset');
hcb=getappdata(r.multihandle,'ColorbarHandle');
if isempty(hcb)
 cbredge=0;
 cbledge=0;
 cbbedge=0;
 cbtedge=0;
 orientation='vertical';
 if nrows==1; 
  cbspan=1.0; 
 end
elseif ishandle(hcb)
 gcbunits=get(hcb,'Units');
 set(hcb,'Units',units);
 cbposition=get(hcb,'Position');
 cbtightinset=ridgepack_cbextent(r.multihandle);
 set(hcb,'Units',gcbunits);
 orientation=getappdata(r.multihandle,'ColorbarOrientation');
 if strcmp(orientation,'vertical')
  if nrows==1; 
   cbspan=1.0; 
  end
  cbredge=minsep+cbposition(3)+cbtightinset(3);
  cbledge=0;
  cbbedge=max(0,(cbtightinset(2)+cbposition(4))-(cbposition(4)/cbspan));
  cbtedge=cbbedge;
 elseif strcmp(orientation,'horizontal')
  if ncols==1; 
   cbspan=1.0; 
  end
  cbredge=max(0,(cbtightinset(1)+cbposition(3))-(cbposition(3)/cbspan));
  cbledge=cbredge;
  cbbedge=minsep+cbposition(4)+cbtightinset(2);
  cbtedge=0;
 else
  error('Incorrect orientiation specification for global colorbar')
 end
else
 error('Accessing colorbar handle')
end

% if title axes has a legend, get required edge gap and set legend position
lph=getappdata(r.titlehandle,'LegendPeerHandle');
tledge=0; rledge=0; bledge=0; lledge=0;
if ~isempty(lph)
 titledata=get(r.titlehandle,'UserData');
 location=get(lph,'Location');
 if ~isfield(titledata,'legendloc')
  titledata.legendloc=[];
 end
 if ~isempty(strfind(location,'outside')) | ~isempty(strfind(location,'none')) 
  lposition=get(lph,'Position');
  if ~isempty(strfind(location,'north')) | strcmp(titledata.legendloc,'north')
   titledata.legendloc='north';
   tledge=lposition(4);
  elseif ~isempty(strfind(location,'east')) | strcmp(titledata.legendloc,'east')
   titledata.legendloc='east';
   rledge=lposition(3);
  elseif ~isempty(strfind(location,'west')) | strcmp(titledata.legendloc,'west')
   titledata.legendloc='west';
   lledge=lposition(3);
  else
   if isempty(strfind(location,'south')) & ~strcmp(titledata.legendloc,'south')
    disp('Defaulting to South Position')
   end
   titledata.legendloc='south';
   bledge=lposition(4);
  end
 elseif debug
  disp('Hmmm... global title axes has an interior legend')
 end
 set(r.titlehandle,'UserData',titledata)
end

% make room for title header if required
titledata=get(r.titlehandle,'UserData');
if ~isempty(titledata) & isfield(titledata,'overshoot')
 tovershoot=titledata.overshoot;
else
 tovershoot=0;
end

% Set new x gap values based on the tight insets
r.gap.x(1)=minedge+cbledge+lledge+max([r.tightinset(:,1,1)]);
if ncols>1
 for i=2:ncols
  r.gap.x(i)=max([r.tightinset(:,i-1,3)+r.tightinset(:,i,1)]);
 end
end
r.gap.x(end)=minedge+cbredge+rledge+max([r.tightinset(:,end,3)]);

% Set new y gap values based on the tight insets
r.gap.y(1)=tovershoot+minedge+cbtedge+tledge+titleinset(4)+max([r.tightinset(1,:,4)]);
if nrows>1
 for i=2:nrows
  r.gap.y(i)=max([r.tightinset(i-1,:,2)+r.tightinset(i,:,4)]);
 end
 if ncols==1
  r.gap.y(2:nrows)=max(r.gap.y(2:nrows));
 end
end
r.gap.y(end)=minedge+cbbedge+bledge+max([r.tightinset(end,:,2)]);

% Adjust extents if there is over/undershoot in x direction
spanx=rx-sum(r.gap.x(:))-sum(r.extent.x(:));
if ncols>1 
 for i=1:length(r.extent.x)
  r.extent.x(i)=r.extent.x(i)+spanx/ncols;
 end
 for i=1:length(r.extent.y)
  r.extent.y(i)=r.extent.y(i)+spanx/ncols;
 end
else
 for i=1:length(r.extent.x)
  r.extent.x(i)=r.extent.x(i)+spanx;
 end 
end

% Adjust extents if there is over/undershoot in y direction
% Note that because correction to the y-extent is carried out
% after the x-extent, the gap above the axes cluster must be
% greater than or equal to the gap below, to ensure correct
% alignment. This is why the colorbar edge is added to the top
% and bottom of the figure in gap calculations above.
spany=ty-sum(r.gap.y(:))-sum(r.extent.y(:));
if spany<0
 for i=1:length(r.extent.y)
  r.extent.y(i)=r.extent.y(i)+spany/nrows;
 end
 if ncols>1
  for i=1:length(r.extent.x)
   r.extent.x(i)=r.extent.x(i)+spany/nrows;
  end
 end
end

% calculate final residual space on the multiplot
spanx=rx-sum(r.gap.x(:))-sum(r.extent.x(:));
spany=ty-sum(r.gap.y(:))-sum(r.extent.y(:));

% determine left top coordinates of axes being plotted
if col>1
 x=spanx/2+sum(r.gap.x(1:col))+sum(r.extent.x(1:col-1))-r.dead.x(row,col);
else
 x=spanx/2+sum(r.gap.x(1:col))-r.dead.x(row,col);
end
y=ty-spany/2-sum(r.gap.y(1:row))-sum(r.extent.y(1:row))-r.dead.y(row,col);

% change position of the axes
if isempty(ha)
 ha=axes('Position',[x/rx y/ty r.extent.x(col)/rx r.extent.y(row)/ty]);
 haunits=get(ha,'Units');
 set(ha,'Units',units);
else
 set(ha,'Units',units);
 set(ha,'Position',[x y r.extent.x(col) r.extent.y(row)]);
end

% adjust axes for maps
NewPosition=ridgepack_axpos(ha);
FullPosition=getpixelposition(ha);

% convert new positions in pixels to normalized extent
r.extent.x(col)=NewPosition(3)/widthfigpixels;
r.extent.y(row)=NewPosition(4)/heightfigpixels;

r.dead.x(row,col)=(NewPosition(1)-FullPosition(1))/widthfigpixels;
r.dead.y(row,col)=(NewPosition(2)-FullPosition(2))/heightfigpixels;

% apply this as the best first guess for all axes 
% not yet created in the same row or column
for ips=1:nrows; for jps=1:ncols;
 if isnan(r.handle(ips,jps)) & (row==ips | col==jps)
  r.extent.x(jps)=r.extent.x(col);
  r.extent.y(ips)=r.extent.y(row);
  r.dead.x(ips,jps)=r.dead.x(row,col);
  r.dead.y(ips,jps)=r.dead.y(row,col);
 end
end; end;

% determine left top coordinates of axes being plotted
if col>1
 x=spanx/2+sum(r.gap.x(1:col))+sum(r.extent.x(1:col-1))-r.dead.x(row,col);
else
 x=spanx/2+sum(r.gap.x(1:col))+r.dead.x(row,col);
end
y=ty-spany/2-sum(r.gap.y(1:row))-sum(r.extent.y(1:row))-r.dead.y(row,col);

% reposition current frame
set(ha,'Position',[x y r.extent.x(col) r.extent.y(row)]);

% position title axes with slight over shoot at top of axes
width=sum(r.gap.x(2:end-1))+sum(r.extent.x(1:ncols));
height=sum(r.gap.y(2:end-1))+sum(r.extent.y(1:nrows));
if col>1
 x=(spanx/2)+r.gap.x(1);
end
y=ty-(spany/2)-r.gap.y(1)-height;

% position title axes
set(r.titlehandle,'Units',units);
set(r.titlehandle,'Position',[x y width height+max(r.tightinset(1,:,4))+tledge])

% reposition the global legend if required and accomodate top legends
if ishandle(lph)
 if strcmp(titledata.legendloc,'north')
   lposition(2)=y+height+max(r.tightinset(1,:,4));
 elseif strcmp(titledata.legendloc,'east')
   lposition(1)=x+width+cbredge+max(r.tightinset(:,end,3));
 elseif strcmp(titledata.legendloc,'west')
   lposition(1)=x-lledge-max(r.tightinset(:,1,1));
 elseif strcmp(titledata.legendloc,'south')
   lposition(2)=y-bledge-cbbedge-max(r.tightinset(end,:,2));
 else
   error('Legend position error')
 end
 set(lph,'Position',lposition,'Box','off');
end

% adjust colorbar axis to a fraction of height or width of multiplot by
% changing the vertical or horizontal extent of the associated axes
if strcmp(orientation,'vertical')
 y=y+(1-cbspan)*height/2;
 height=cbspan*height;
 width=width+max(r.tightinset(:,end,3));
elseif strcmp(orientation,'horizontal')
 x=x+(1-cbspan)*width/2;
 width=cbspan*width;
 y=y-max(r.tightinset(end,:,2));
else
 error('orientiation specification error for global axes')
end
set(r.multihandle,'Position',[x y width height]);

% add notation to the multifigure structure or replot if needed
if nargin>=7 & ischar(notation)
 r.notation{row,col}=notation;
end
if ~isempty(char(r.notation{row,col})) & ishandle(r.notationhandle(row,col))
 delete(r.notationhandle(row,col))
end
if ~isempty(char(r.notation{row,col})) & plotnotation
 NewPosition=ridgepack_axpos(ha);
 posratio=NewPosition(4)/NewPosition(3);
 if fontbox==5 % left 1/3 down left side
  r.notationhandle(row,col)=text(min(0.02,0.05*posratio),0.75,char(r.notation{row,col}),...
                                'HorizontalAlignment','left',...
                                'VerticalAlignment','middle',...
                                'Units','normalized',...
				'FontName','Helvetica',...
                                'Interpreter','Tex',...
                                'FontSize',fontsize,'Color',fontcol,...
                                'Parent',ha,'Interpreter','Tex');
 elseif fontbox==4 % upper right
  r.notationhandle(row,col)=text(min(0.98,0.96*posratio),0.97,char(r.notation{row,col}),...
                                'HorizontalAlignment','right',...
                                'VerticalAlignment','top',...
                                'Units','normalized',...
				'FontName','Helvetica',...
                                'Interpreter','Tex',...
                                'FontSize',fontsize,'Color',fontcol,...
                                'Parent',ha,'Interpreter','Tex');
 elseif fontbox==3 % lower left
  r.notationhandle(row,col)=text(min(0.02,0.05*posratio),0.03,char(r.notation{row,col}),...
                                'HorizontalAlignment','left',...
                                'VerticalAlignment','bottom',...
                                'Units','normalized',...
				'FontName','Helvetica',...
                                'Interpreter','Tex',...
                                'FontSize',fontsize,'Color',fontcol,...
                                'Parent',ha,'Interpreter','Tex');
 elseif fontbox==2 % lower right
  r.notationhandle(row,col)=text(min(0.98,0.96*posratio),0.03,char(r.notation{row,col}),...
                                'HorizontalAlignment','right',...
                                'VerticalAlignment','bottom',...
                                'Units','normalized',...
				'FontName','Helvetica',...
                                'Interpreter','Tex',...
                                'FontSize',fontsize,'Color',fontcol,...
                                'Parent',ha,'Interpreter','Tex');
 elseif fontbox==1 % upper right with box around it
  r.notationhandle(row,col)=text(max(0.02,0.05*posratio),0.97,char(r.notation{row,col}),...
                                'HorizontalAlignment','left',...
                                'VerticalAlignment','top',...
                                'Units','normalized',...
				'FontName','Helvetica',...
                                'Interpreter','Tex',...
                                'FontSize',fontsize,'Color',fontcol,...
				'BackgroundColor',[1 1 1],'EdgeColor','k',...
                                'Parent',ha,'Interpreter','Tex');
 else % upper right
  r.notationhandle(row,col)=text(max(0.02,0.05*posratio),0.97,char(r.notation{row,col}),...
                                'HorizontalAlignment','left',...
                                'VerticalAlignment','top',...
                                'Units','normalized',...
				'FontName','Helvetica',...
                                'Interpreter','Tex',...
                                'FontSize',fontsize,'Color',fontcol,...
                                'Parent',ha,'Interpreter','Tex');
 end
end

% reset units
set(r.titlehandle,'Units',tunits);
set(r.multihandle,'Units',gunits);
set(ha,'Units',haunits);

% set axis font size multiplier
if (nrows>6 | ncols>9)
 set(ha,'LabelFontSizeMultiplier',0.60)
 set(ha,'TitleFontSizeMultiplier',0.60)
elseif (nrows>3 | ncols>5)
 set(ha,'LabelFontSizeMultiplier',0.80)
 set(ha,'TitleFontSizeMultiplier',0.80)
end

% set axis handle if it is empty
r.handle(row,col)=ha;

% update the appdata for the multiplot
setappdata(hf,'MultiplotConfiguration',r)

% remove axis handle if no output requested
if nargout==0; clear ha; end

if debug; disp(['...Leaving ',mfilename]); end

