function [hcb]=ridgepack_cbpos(h,orientation)

% ridgepack_cbpos - Colorbar positioning function for native icepack colorbars
%
% function [hcb]=ridgepack_cbpos(h,orientation)
% 
% This function positions the colorbar axes and the main axes with which
% it is associated. It has two inputs, and the output is the colorbar 
% axes handle.
%
% Inputs:
% h	      - handle of axes for which a colorbar is required.
% orientation - specified as either 'vertical' colorbar on the right
%               side of the main axes, or a 'horizontal' colorbar 
%               underneath the main axes.
%
% Output:
% hcb	      - colorbar handle.
%
% Andrew Roberts, Naval Postgraduate School, March 2018  (afrobert@nps.edu)

global debug;
if debug; disp(['Entering ',mfilename,'...']); end

% Colorbar positioning and axis limit constants
cbaspect=0.025; % default aspect ratio (height) of vertical (horizontal) bar
cmap=colormap; % gives the color limits of the axis
minsep=0.01; % minumum separation of colorbar from axes in normalized units
             % (must be same as in ridgepack_multipos.m)
minwidth=0.015; % minimum width of colorbar in normalized units

% Find axis and colorbar handles
if nargin==0
 error('No axis handle input')
elseif ishandle(h) && any(strcmp(get(h,'Type'),{'axes'}))
 axhandle=get(h);
 hcb=getappdata(h,'ColorbarHandle');
 hdl=getappdata(h,'ColorbarListener');
else
 error('Input is not an axes handle')
end


% check inputs
if isempty(hcb) && (nargin<2 || ~ischar(orientation))
 error('Colorbar orientation missing')
elseif isempty(hcb)
 if debug; disp('Generating new Colorbar'); end
elseif ishandle(hcb) 
 orientation=getappdata(h,'ColorbarOrientation');
elseif ischar(hcb) && strcmp(hcb,'ColorbarDeleted')
 % reached here if colorbar has been deleted
 return
elseif ~ishandle(hcb)
 ridgepack_cbdelete(h);
 return
else
 error('Axis does not contain colorbar axis handle in expected format')
end

% Check orientation specification
if ~any(strcmp(orientation,{'vertical','horizontal'}))
 error('Colorbar orientation is incorrect')
end

% change aspect ratio and minwidth if colorbar spans whole figure
axpos=get(h,'Position'); % get parent axes positions
if strcmp(orientation,'vertical') & axpos(4)>0.70
 cbaspect=0.019; 
 minwidth=0.012; 
elseif strcmp(orientation,'horizontal') & axpos(3)>0.70
 cbaspect=0.018; 
 minwidth=0.012; 
end

% get extent of figure in pixels
hfpos=getpixelposition(get(h,'Parent'));
widthfigpixels=hfpos(3);
heightfigpixels=hfpos(4);

% Find Position in pixels of main axes and get the tight inset
Position=getpixelposition(h);

% Calculation VSeparation and HSeparation parameters
% This is further adjusted by adding the axes tightinset later
VSeparation=minsep*widthfigpixels;
HSeparation=minsep*heightfigpixels;

% Reposition axes if no colorbar exists
if isempty(hcb) 
 if strcmp(orientation,'vertical')
  Position(3)=Position(3)-VSeparation-Position(4)*cbaspect;
 elseif strcmp(orientation,'horizontal')
  Position(4)=Position(4)-HSeparation-Position(3)*cbaspect;
  Position(2)=Position(2)+HSeparation+Position(3)*cbaspect;
 else
  error('Incorrect orientiation specification')
 end
end

% set position regardless of whether it has changed
setpixelposition(h,Position)

% Real position of axes (different from 'position' in the axes handle)
NewPosition=ridgepack_axpos(h);

% Generate colorbar position 
if strcmp(orientation,'vertical')

 CBPosition(3)=NewPosition(4)*cbaspect;
 CBPosition(1)=NewPosition(1)+NewPosition(3)+VSeparation;
 CBPosition(4)=NewPosition(4);
 CBPosition(2)=NewPosition(2);

 xlim=[-1 2]; 
 ylim=0.5+[0 length(cmap)];

elseif strcmp(orientation,'horizontal')

 CBPosition(3)=NewPosition(3);
 CBPosition(1)=NewPosition(1);
 CBPosition(4)=NewPosition(3)*cbaspect;
 CBPosition(2)=NewPosition(2)-CBPosition(4)-HSeparation;

 xlim=0.5+[0 length(cmap)];
 ylim=[-1 2]; 

else

 error('Incorrect orientiation specification')

end

% Create new colorbar or reposition existing colorbar
if isempty(hcb)

 % Generate colorbar axes and listener
 hcb=axes('xlim',xlim,'ylim',ylim,'box','on','Layer','top');
 setpixelposition(hcb,CBPosition)
 set(hcb,'Xtick',[],'XTickLabel',[],'Ytick',[],'YTickLabel',[])
 set(hcb,'tag','Colorbar','ActivePositionProperty','Position')

 % Move main figure axes to accommodate colorbar
 setpixelposition(h,Position)
 set(h,'ActivePositionProperty','Position')
 set(h,'Tag','MainAxes')

 %h = myclass();

 % Create axes listener to automatically reposition the colorbar with figure position
 hdl=addlistener(h,{'XLim','YLim','ZLim','Position','DataAspectRatio',...
                    'PlotBoxAspectRatio','DataAspectRatioMode',...
                    'PlotBoxAspectRatioMode',...
                    'OuterPosition','XLabel','YLabel','ZLabel'},'PostSet',@ridgepack_cbfix);

%                    'PlotBoxAspectRatioMode','YTickLabel','XTickLabel','ZTickLabel',...
%                    'OuterPosition','XLabel','YLabel','ZLabel'},'PostSet',@ridgepack_cbfix);

 % Add colorbar axes, listener and orientation to main axes  app data
 setappdata(h,'ColorbarHandle',hcb)
 setappdata(h,'ColorbarListener',hdl)
 setappdata(h,'ColorbarOrientation',orientation)

 % Set delete function remove colorbar if the main axes are deleted
 set(h,'DeleteFcn','ridgepack_cbdelete(gcbo)')

elseif ishandle(hcb)

 % resposition colorbar 
 setpixelposition(hcb,CBPosition)
 set(hcb,'box','on')

else

 error('No colorbar positioning or repositioning occurred')

end

% Get position information
CBPosition=get(hcb,'Position');

% Set fontsize
if strcmp(orientation,'vertical')
 fontsize=min(max(8,12*CBPosition(4).^(1/3)),10);
elseif strcmp(orientation,'horizontal')
 fontsize=min(max(8,12*CBPosition(3).^(1/3)),10);
end
set(hcb,'FontSize',fontsize);

% Set fontsize of children
child=get(hcb,'Children');
if ~isempty(child)
 for i=1:length(child)
  if ishandle(child(i)) && any(strcmp(get(child(i),'Type'),{'text'}))
   set(child(i),'FontSize',fontsize);
  end
 end
end

% Get position information
set(h,'ActivePositionProperty','Position')
AXPosition=get(h,'Position');
AXTightInset=get(h,'TightInset');
CBPosition=get(hcb,'Position');
CBTightInset=ridgepack_cbextent(h);

% Extract information on sharing a colorbar
share=getappdata(h,'ColorbarShare');

% Adjust colorbar width using normalized position and tightinset
if strcmp(orientation,'vertical')
 CBPosition(3)=max(minwidth,cbaspect*CBPosition(4));
 if AXTightInset(3)==0 & strcmpi(get(h,'YAxisLocation'),'right') & ...
    (~isempty(get(get(h,'YLabel'),'String')) & ~isempty(get(h,'YTickLabel')))
  CBPosition(1)=min(CBPosition(1)+0.075,1-CBPosition(3)-CBTightInset(3));
 elseif AXTightInset(3)==0 & strcmpi(get(h,'YAxisLocation'),'right') & ...
    (~isempty(get(get(h,'YLabel'),'String')) | ~isempty(get(h,'YTickLabel')))
  CBPosition(1)=min(CBPosition(1)+0.04,1-CBPosition(3)-CBTightInset(3));
 elseif AXTightInset(3)==0 & ~isempty(get(h,'XTickLabel')) 
  CBPosition(1)=min(CBPosition(1)+AXTightInset(3),1-CBPosition(3)-CBTightInset(3));
 else
  CBPosition(1)=min(CBPosition(1)+AXTightInset(3),1-CBPosition(3)-CBTightInset(3));
 end
 if ~isempty(share) 
  CBPosition(4)=1.75*CBPosition(4);
  CBPosition(2)=CBPosition(2)-CBPosition(4)/2-(AXTightInset(4)/2)-AXTightInset(2);
 end
elseif strcmp(orientation,'horizontal')
 CBPosition(4)=max(minwidth,cbaspect*CBPosition(3)*widthfigpixels/heightfigpixels);
 if AXTightInset(2)==0 & strcmpi(get(h,'XAxisLocation'),'bottom') & ...
    (~isempty(get(get(h,'XLabel'),'String')) & ~isempty(get(h,'XTickLabel')))
  CBPosition(2)=max(0,CBPosition(2)-0.075);
 elseif AXTightInset(2)==0 & strcmpi(get(h,'XAxisLocation'),'bottom') & ...
    (~isempty(get(get(h,'XLabel'),'String')) | ~isempty(get(h,'XTickLabel')))
  CBPosition(2)=max(0,CBPosition(2)-0.04);
 else
  CBPosition(2)=max(0,CBPosition(2)-AXTightInset(2)-minsep);
 end
 if ~isempty(share) 
  %CBPosition(1)=CBPosition(1)+0.25*CBPosition(3)/2+(AXTightInset(1)/2)+AXTightInset(3); 
  CBPosition(1)=CBPosition(1)+0.25*CBPosition(3)/2+AXTightInset(1)+AXTightInset(3); 
  CBPosition(3)=1.75*CBPosition(3);
 end
end

% Set colorbar axes position
set(hcb,'Position',CBPosition);

if debug; disp(['...Leaving ',mfilename]); end

