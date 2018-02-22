function ridgepack_clearax(axisxy,keeptick)

% ridgepack_clearax - Clear the current axes labels and tick marks
%
% function ridgepack_clearax(axisxy,keeptick)
%
% This function clears tick marks and axis labels from the
% current axes and centers the axes positions accordingly.
%
% Inputs:
% axisxy   - Specifies if only the x or y axes are to be cleared
%            This is a character variables with allowable values
%            of 'x' or 'y'.  If it is omitted, then all axes are 
%            cleared.
% keeptick - Optional argument, keep the tick marks when set to 1. 
%            The default is zero.
% 
% Ridgepack Version 1.0
% Andrew Roberts, Naval Postgraduate School, March 2018 (afrobert@nps.edu)
%

global debug;
if debug; disp(['Entering ',mfilename,'...']); end


if nargin==0
 clearx=1;
 cleary=1;
elseif ~ischar(axisxy)
 error('Input must be a string: ''x'' or ''y''')
elseif strcmp(axisxy,'x')
 clearx=1;
 cleary=0;
elseif strcmp(axisxy,'y')
 clearx=0;
 cleary=1;
elseif strcmp(axisxy,'both')
 clearx=1;
 cleary=1;
else
 clearx=0;
 cleary=0;
end

if nargin<2
 keeptick=0;
elseif nargin==3 & (~isnumeric(keeptick) | (keeptick<0 | keeptick>1))
 error('keeptick should be either 0 or 1') 
end

if clearx
 xlabel('');
 if keeptick
  set(gca,'XTickLabel',[])
 else
  set(gca,'Xtick',[],'XTickLabel',[])
 end
end

if cleary
 ylabel('');
 if keeptick 
  set(gca,'YTickLabel',[])
 else
  set(gca,'Ytick',[],'YTickLabel',[])
 end
end

drawnow;

if debug; disp(['...Leaving ',mfilename]); end


