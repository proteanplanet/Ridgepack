function ridgepack_multialign(hf,titletext,notatsize,notatcol,notatbox)

% ridgepack_multialign - Align all axes on a figure constructed using ridgepack_multiplot
%
% function ridgepack_multialign(hf,titletext,notatsize,notatcol,notatbox)
%
% This function aligns axes on a figure constructed using ridgepack_multiplot. It 
% is called automatically when there is a change of axes that are constructed 
% with ridgepack_multiplot, and must be called manually to finalize the figure. 
% In order to use this function, you must first create axes using the 
% ridgepack_multiplot function. If complete alignment is not achieved when this is 
% first executed, it may be executed or as many iterations as desired, although
% a need to run ridgepack_multialign multiple times to line up axes is rare.
%
% INPUT:
%
% hf        - Figure handle (if omitted, the current figure is adopted)
% titletext - Title text of multiple axes figure (optional).  A title added 
%             can later be removed by including titletext as an empty string 
%             (i.e. set titletext='') by re-running ridgepack_multialign.
% notatsize - Font size of the text used for notations (typically 8-14; optional).
% notatcol  - Font color of the text used for notations (RBG vector; optional).
% notatbox  - If set to 1, provides a box around the notation, if set to 2,
%             puts notation in lower right, if set to 3, puts it in the lower
%             left, and if set to 4, puts it in the upper right. If set to 5
%             it places the number 1/3 down the left side.
%
% Example:
% ridgepack_multiplot(1,2,1,1,'a') , ridgepack_polarm('seaice');
% ridgepack_multiplot(1,2,1,2,'b') , ridgepack_polarm('seaice');
% ridgepack_multialign(gcf,'This is the title')
% Once all calls have been made to ridgepack_multiplot, a final call of ridgepack_multialign 
% brings them into alignement and uses all available space on the figure. 
%
% Ridgepack Version 1.0
% Andrew Roberts, Naval Postgraduate School, March 2018 (afrobert@nps.edu)


global debug;
if debug; disp(['Entering ',mfilename,'...']); end

if nargin<1; 
 hf=gcf; 
elseif ~strcmpi(get(hf,'type'),'figure')
 error('handle provided is not a figure handle')
end

% check font size and font color
if nargin>2 & ~isnumeric(notatsize)
 error('Font size for notation must be a number')
elseif nargin>3 & ~isnumeric(notatcol)
 error('Font color for notation must be a number')
elseif nargin>3 & length(notatcol)~=3
 error('Font color must be a three element RGB vector')
elseif nargin<=2
 notatsize=12;
 notatcol=[0 0 0];
elseif nargin<=3
 notatcol=[0 0 0];
end

if nargin>=5 & ~isnumeric(notatbox)
 error('notatbox is not a number')
elseif nargin<5
 notatbox=0;
end

if ishandle(hf)

 r=getappdata(hf,'MultiplotConfiguration');

 if nargin>=2 & ischar(titletext) & ~isempty(r)


  % assign global title if requested, and determine overshoot over top of 
  % figure so as to make space at top of figure if required
  if ishandle(r.titlehandle)

   titledata=get(r.titlehandle,'UserData');
   
   ha=gca;
   set(hf,'CurrentAxes',r.titlehandle)
   title(titletext,'FontWeight','normal')
   titleaxpos=get(r.titlehandle,'Position');
   titleextent=get(r.titlehandle,'TightInset');

   titledata.type='title';
   if ~isfield(titledata,'overshoot')
    titledata.overshoot=max(0,titleaxpos(2)+titleaxpos(4)+titleextent(4)-1);
   elseif nargin>=2 & strcmp(titletext,'')
    titledata=rmfield(titledata,'overshoot');
   end

   set(r.titlehandle,'UserData',titledata)
   set(hf,'CurrentAxes',ha)

  end

 elseif nargin>=2 & ~ischar(titletext) & ~isempty(r)

  error('Cannot recognize title input')

 end

 if ~isempty(r) & isstruct(r) & ...
     isfield(r,'handle') & isfield(r,'extent') & ...
     isfield(r.extent,'x') & isfield(r.extent,'y')

  nrows=length(r.extent.y);
  ncols=length(r.extent.x);

  % now commence alignement 
  oldx=0; oldy=0;
  pos=get(r.multihandle,'Position');
  newx=pos(1); newy=pos(2);
  k=0;
  postolerance=0.000025;

  while (abs(newx-oldx)>postolerance | abs(newy-oldy)>postolerance | k<3) & k<10

   oldx=newx; oldy=newy;

   for row=1:nrows
   for col=1:ncols
     hac=r.handle(row,col);
     if ishandle(hac) 
      ridgepack_multipos(hf,hac,nrows,ncols,row,col)
     end
   end
   end

   % get configuration
   r=getappdata(hf,'MultiplotConfiguration');
   pos=get(r.multihandle,'Position');
   newx=pos(1); newy=pos(2);
   k=k+1;

  end

  % final alignment, adding notation text (not added until now
  % because printed graphics will include each time the 
  % graphics are plotted, even if the user interface does not.
  for row=1:nrows
  for col=1:ncols
    hac=r.handle(row,col);
    if ishandle(hac) 
     ridgepack_multipos(hf,hac,nrows,ncols,row,col,1,notatsize,notatcol,notatbox)
    end
  end
  end


  if k==10; 
   disp('Reached ten iterations of multiplot alignment')
  end

 else

  error('Unable to access multiplot information')

 end

else

 error('Figure handle error')

end

drawnow

if debug; disp(['...Leaving ',mfilename]); end


