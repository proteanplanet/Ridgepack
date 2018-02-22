function ridgepack_vecref(h,scalelength,arrow,refval,units,veccol,vecwidth,accenture,lx,ly)

% ridgepack_vecref - Creates reference vector for a 2D vector diagram
%
% function ridgepack_vecref(h,scalelength,arrow,refval,units,veccol,vecwidth,accenture)
%
% This function creates a reference vector to accompany a (non color-coded)
% vector field plot. As part of the vector key creation, a data structure
% that describes information required to create the key is added to the 
% app data of the main axes to which the key is assigned.
%
% Inputs:
% h	      - main axes handle on which the vector key is to be added
% scalelength - the scalelength of all vectors, such that every vector on 
%               the figure is set by [u,v]*scalelength/refval
% arrow       - describes the shape of the arrow to be used in the reference
%               vector key (should be the same as in main figure)
% refval      - reference value represented in the key
% units	      - units to be plotted as part of the reference vector key
% veccol      - color of the reference vector key
% vecwidth    - line width of vector to be plotted
% accenture   - shape of vectors (1=standard arrow)
% lx,ly       - vector locations from ridgepack_quiverref
%
% Ridgepack Version 1.0
% Andrew Roberts, Naval Postgraduate School, March 2018 (afrobert@nps.edu)

global debug;
if debug; disp(['Entering ',mfilename,'...']); end

% Remove reference vector app data if needed
rcc=getappdata(h,'ReferenceVectorStructure');
dcc=getappdata(h,'ReferenceVectorDeleted');
if ~isempty(rcc)
 rmappdata(h,'ReferenceVectorStructure');
elseif dcc
 if debug; disp('Reference vector key deleted'); end
 setappdata(h,'ReferenceVectorDeleted',true);
 return
else
 if debug; disp('Creating vector reference key'); end
end

% Ensure that the data aspect ratio is 1:1 in the x:y direction
ratio=get(h,'DataAspectRatio');
if ~all(ratio==[1 1 1])
 disp('Data Aspect Ratio must be 1:1')
 set(h,'DataAspectRatio',[1 1 1]);
end

% Get the reference text string, formatted to powers of ten if required
if refval < 0.01 | refval > 100
 factor=floor(log10(refval));
 reftext=[num2str(refval/(10^factor)),' \times 10^{',num2str(factor),'} ',units,' '];
else
 reftext=[num2str(refval),' ',units,' '];
end

% Get the current axis limits
xlimp=xlim(h); xp1=xlimp(1); xp2=xlimp(2);
ylimp=ylim(h); yp1=ylimp(1); yp2=ylimp(2);

% calculate length of 2D ticks in local coordinates and add
% this as a buffer, so that it is consistent with MATLAB legends
% Note that adding 0.0025 "seems to work" when position(4)>=position(3)
% when printing out the figures
if ismap(h)
 tl=0; % THIS SETS THE VECTOR BOX AGAINST THE FRAME EDGE
else
 gposition=get(h,'Position');
 position=ridgepack_axpos(h);
 ticklength=get(h,'TickLength');
 if position(3)>position(4)
  ticklength(1)=ticklength(1);
  tl=ticklength(1)*diff(xlimp)./max(gposition(3:4));
 elseif position(4)>=position(3)
  ticklength(1)=(ticklength(1)+0.0025);
  tl=ticklength(1)*diff(ylimp)./max(gposition(3:4));
 end
end

% set padding around the reference vector
padx=diff(xlimp)/50;
pady=diff(ylimp)/100;

% Set x position of reference vector
xend=xp2-padx-tl;
xstart=xend-scalelength;

% adjust fontsize for vertical height of plot
%figout=get(h,'OuterPosition');
%fontsize=min(11,max(8,10*sqrt((figout(3).^2)+(figout(4).^2))/sqrt(2)));
fontsize=get(gca,'fontsize');

% Plot reference text in lower right hand corner to get dimensions
ht=text(xstart-padx,yp1+pady+tl,reftext,'Visible','off','Parent',h,...
        'FontSize',fontsize,'VerticalAlignment','Middle',...
        'HorizontalAlignment','Right','Interpreter','tex',...
        'Clipping','off');
textextent=get(ht,'Extent');
delete(ht);

% Find place of lowest density of vectors to place key
bxr=xp2-tl;
bxl=bxr-scalelength-textextent(3)-2*padx;
byb=yp1+tl;
byt=byb+textextent(4)+2*pady;

txr=xp2-tl;
txl=txr-scalelength-textextent(3)-2*padx;
tyt=yp2-tl;
tyb=tyt-textextent(4)-2*pady;

lx=lx(~isnan(lx));
ly=ly(~isnan(ly));
bidx=find(lx>bxl & lx<bxr & ly>byb & ly<byt);
tidx=find(lx>txl & lx<txr & ly>tyb & ly<tyt);

% Draw patch over area of vector key 
z=0.001;
if length(bidx)>length(tidx) 
 hp=patch([txl; txl; txr; txr],[tyb; tyt; tyt; tyb],[z; z; z; z],'w','Parent',h);
 yref=(tyb+tyt)/2;
else
 hp=patch([bxl; bxl; bxr; bxr],[byb; byt; byt; byb],[z; z; z; z],'w','Parent',h);
 yref=(byb+byt)/2;
end
uistack(hp,'top');

% Redraw reference text on top of patch
ht=text(xstart,yref,2.1,reftext,'Parent',h,'FontSize',fontsize,...
         'VerticalAlignment','Middle','HorizontalAlignment','Right',...
         'Interpreter','tex');

% Set y position of reference vector
%yend=yb+(textextent(4)+2*pady)/2;
%ystart=yend;
yend=yref;
ystart=yref;

% Get x coordinates of reference vector plotted
lx = [xstart; ...
      xstart+(1-arrow*accenture)*(xend-xstart); ...
      xend-arrow*scalelength; ...
      xend; ...
      xend-arrow*scalelength; ...
      xstart+(1-arrow*accenture)*(xend-xstart); ...
      NaN];


% Get y coordinates of reference vector plotted
ly = [ystart; ...
      ystart+(1-arrow*accenture)*(yend-ystart); ...
      yend+arrow*(arrow*scalelength); ...
      yend; ...
      yend-arrow*(arrow*scalelength); ...
      ystart+(1-arrow*accenture)*(yend-ystart); ...
      NaN];

% Get z coordinates of reference vector
lz = 2*ones(size(ly));

% Plot the reference vector
hl=line(lx,ly,lz,'Color',veccol,'Parent',h,'LineWidth',vecwidth);

% create structure to attach to axes data
rc.scalelength=scalelength;
rc.arrow=arrow;
rc.refval=refval;
rc.units=units;
rc.veccol=veccol;
rc.vecwidth=vecwidth;
rc.accenture=accenture;
rc.lx=lx;
rc.ly=ly;
rc.ht=ht;
rc.hp=hp;
rc.hl=hl;

% Add structure as appdata to current axes
setappdata(h,'ReferenceVectorStructure',rc)
setappdata(h,'ReferenceVectorDeleted',false);

% Create a new listener for this key, write to app data
hdl=addlistener(h,{'Position','DataAspectRatio','PlotBoxAspectRatio',...
                   'DataAspectRatioMode','PlotBoxAspectRatioMode',...
                   'XLim','YLim','ZLim','Visible'},'PostSet',@ridgepack_vecfix);
setappdata(h,'ReferenceVectorListener',hdl);

if debug; disp(['...Leaving ',mfilename]); end

