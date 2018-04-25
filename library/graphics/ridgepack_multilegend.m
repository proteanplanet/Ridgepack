function [lh]=ridgepack_multilegend(objecthandles,legendtext,location)

% ridgepack_multilegend - Create a global legend for a figure created with ridgepack_multiplot
%
% function [lh]=ridgepack_multilegend(objecthandles,legendtext,location)
%
% This function creates a legend for a figure created with ridgepack_multiplot, where
% a large legend may be required to improve legibility, or where it is not 
% apparent that a legend in a single set of axes applies to axes plotted in the
% figure. The legend will always be plotted outside of the cluster of axes 
% created with ridgepack_multiplot, and the North and South options produce horizontal
% legends, while the East and West legends are vertical.  The default is to
% not include a box around the legend, but this can be changed by using the 
% legend output handle from this function and setting: legend(lh,'boxoff').
%
% INPUT:
%
% objecthandles - Handles created from graphics objects (e.g. plot, contour, etc.)
%                 that you would like represented in the legend.
% legendtext    - A cell array of the text to accompany each legend entry,
%                 corresponding to each object handle provided.
% location      - Desired location of the legend outside the plot region
%                 This can take the options 'North','West','East' and 'South'
%                 and defaults to South.
%
%
% OUTPUT:
%
% lh - this is the handle of the MATLAB legend created. 
%
% Example 1: Adding legends to a multiple graph axes
%
%    clear; clf
%    ridgepack_multiplot(3,1,1,1,'a');             
%    h(1)=plot([5:25],[5:25].^2);
%    ha(1)=gca;
%    ridgepack_clearax('x',1); title('');
%    ridgepack_multiplot(3,1,2,1,'b');
%    h(2)=plot([1:20],sin([1:20]/3),'r.');
%    ha(2)=gca;
%    ridgepack_clearax('x',1); title('');
%    ridgepack_multiplot(3,1,3,1,'c');
%    h(3)=plot([1:20],cos([2:21]/2),'g--');
%    ha(3)=gca;
%    ridgepack_clearax('none',1); title('');
%    linkaxes(ha,'x')
%    ridgepack_multialign(gcf,'Title of plot')
%    ridgepack_multilegend(h,{'test1','test and more','the ultimate test'})
%
% Example 2: Creating a figure with image, map, colorbar, and legend
%
%    clf
%    clear hs
%    nc=ridgepack_clone('r30RB1FIY.cpl.ha.a2xavg_Sa_pslv.1992.nc',{'a2xavg_Sa_pslv'},1)
%    ridgepack_multiplot(1,2,1,1,'a'), ridgepack_image(nc,'a2xavg_nx','a2xavg_ny','a2xavg_Sa_pslv',...
%                                      {},{},[],'linear',0,'horizontal');
%    title('title goes here'); 
%    ridgepack_multicb(gca);
%    ridgepack_multiplot(1,2,1,2,'b'), ridgepack_polarm('seaice')
%    title('measurement sites')
%    hs(1)=plotm(70,180,'r.');
%    hs(2)=plotm(80,180,'g.');
%    ridgepack_multilegend(hs,{'entry 1','entry 2'},'North')
%    ridgepack_multialign(gcf,'This is the title')
%
% Ridgepack Version 1.0
% Andrew Roberts, Naval Postgraduate School, March 2018 (afrobert@nps.edu)


global debug;
if debug; disp(['Entering ',mfilename,'...']); end


if nargin<2
 error('must provide two arguments')
elseif ~iscell(legendtext)
 error('legend text should be a cell array')
elseif length(objecthandles)~=length(legendtext)
 length(objecthandles)
 length(legendtext)
 error('handles and legend text not of equal length')
else
 for i=1:length(objecthandles)
  if ~ishandle(objecthandles(i))
   error('An objecthandle provided is not actually a handle')
  end
 end
end

% font size
if length(legendtext)>5
 fonts=5;
elseif length(legendtext)>4
 fonts=6;
elseif length(legendtext)>3
 fonts=7;
else
 fonts=8;
end


if nargin==2
 location='South';
elseif ~any(strcmpi(location,{'South','North','East','West'}))
 disp('Location error: Defaulting to South location')
 location='South';
end

r=getappdata(gcf,'MultiplotConfiguration');
if isfield(r,'titlehandle') & ishandle(r.titlehandle)
 if strcmpi(location,'North')
  lh=legend(r.titlehandle,objecthandles,legendtext,'FontSize',fonts,...
         'Location','NorthOutside','Orientation','Horizontal');
 elseif strcmpi(location,'West')
  lh=legend(r.titlehandle,objecthandles,legendtext,'FontSize',fonts,...
         'Location','WestOutside','Orientation','Vertical');
 elseif strcmpi(location,'East')
  lh=legend(r.titlehandle,objecthandles,legendtext,'FontSize',fonts,...
         'Location','EastOutside','Orientation','Vertical');
 elseif strcmpi(location,'South')
  lh=legend(r.titlehandle,objecthandles,legendtext,'FontSize',fonts,...
         'Location','SouthOutside','Orientation','Horizontal');
 else
  error('location error')
 end
else
 error('Unable to locate handle for global title axes of multiplot')
end

set(lh,'box','off');

ridgepack_multialign(gcf);

if nargout==0; clear lh; end

drawnow

if debug; disp(['...Leaving ',mfilename]); end


