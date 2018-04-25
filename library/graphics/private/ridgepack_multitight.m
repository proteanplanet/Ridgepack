function [r]=ridgepack_multitight(ha,r,units,mingap,row,col)

% ridgepack_multitight - Gives tight inset margins around axes for ridgepack_multiplot figures
%
% function [r]=ridgepack_multitight(ha,r,units,mingap,row,col)
%
% This functions determines the tight inset margins around axes, required by 
% the function ridgepack_multipos to generate multiple axes figures.  This provides
% different information from the vanilla MATLAB function get(ga,'TightInset')
% because ridgepack_multiplot also needs to know if space should be allowed for
% a colorbar, either vertical or horizontal, surrounding axes. The colorbar
% needs to have been generated using ridgepack_colorbar, not the vanilla MATLAB
% colorbar function.  This function is called by ridgepack_multipos, which in turn
% is called by ridgepack_multiplot and ridgepack_multialign. 
%
% INPUT:
%
% ha      - handle of axis for which the adjusted tight inset is required
% r       - configuration structure for a multiple axes plot (see ridgepack_multipos)
% units   - axes units to be used (e.g. normalized, pixels etc)
% mingap  - minimum allowable gap between different axes
% row,col - row and column of axes for which the tight inset is calculated
%
%
% OUTPUT:
%
% r       - adjusted configuration structure containing adjusted tight insets
%
% Ridgepack Version 1.0
% Andrew Roberts, Naval Postgraduate School, March 2018 (afrobert@nps.edu)


global debug;
if debug; disp(['Entering ',mfilename,'...']); end

if nargin<6
 error('Incorrect inputs for ridgepack_multitight')
end

if ~isempty(ha) & ~isempty(r)

 % get the TightInset information surrounding the axes
 haunits=get(ha,'Units');
 set(ha,'Units',units)
 position=get(ha,'Position');
 set(ha,'Units',haunits);
 tightinset=get(ha,'TightInset');

 % tightinset is unreliable for graphs laid upon each other, so reset to a
 % small value where there is no xlabel and tick marks, seemingly only a problem
 % when a single column of graphs as stacked upone each other.
 if isempty(get(get(ha,'XLabel'),'String')) & isempty(get(ha,'XTickLabel')) & ...
    strcmpi(get(ha,'XAxisLocation'),'bottom')
  tightinset(2)=mingap*2;
 end

 % if a plotyy graph has been used, also use its tightinset data
 hayy=getappdata(ha,'graphicsPlotyyPeer');
 if ~isempty(hayy)
  tightinsetyy=get(hayy,'TightInset');
  tightinset=max(tightinset,tightinsetyy);
 end

 % assign tightinset to multiplot structure
 tightinset(1)=max(mingap,tightinset(1));
 tightinset(2)=max(mingap,tightinset(2));
 tightinset(3)=max(mingap,tightinset(3));
 tightinset(4)=max(mingap,tightinset(4));

 % create extra room for titles between plots
 if row>1 & ~isempty(get(get(ha,'title'),'String')); 
  tightinset(4)=1.5*tightinset(4); 
 end

 %assign initial tightinset
 r.tightinset(row,col,1:4)=tightinset(1:4);

 % determine space required for each individual axes's icepack colorbar
 hahcb=getappdata(ha,'ColorbarHandle');
 haorientation=getappdata(ha,'ColorbarOrientation');
 if ~isempty(hahcb) & ishandle(hahcb) & ~isempty(haorientation)
  share=getappdata(ha,'ColorbarShare');
  hccbunits=get(hahcb,'Units');
  set(hahcb,'Units',units);
  hacbpos=get(hahcb,'Position');
  hacbinset=ridgepack_cbextent(ha);
  set(hahcb,'Units',hccbunits);
  if strcmp(haorientation,'vertical')
   r.tightinset(row,col,3)=hacbpos(1)-position(1)-position(3)+...
                           hacbpos(3)+hacbinset(3)+hacbinset(1);
   if isempty(share) % only account for vertical spacing if not shared
    r.tightinset(row,col,2)=r.tightinset(row,col,2)+hacbinset(2);
    r.tightinset(row,col,4)=r.tightinset(row,col,4)+hacbinset(4);
   end
  elseif strcmp(haorientation,'horizontal')
   r.tightinset(row,col,1)=r.tightinset(row,col,1)+hacbinset(1);
   r.tightinset(row,col,2)=position(2)-hacbpos(2)+hacbinset(2)+hacbinset(4);
   if isempty(share) % only account for vertical spacing if not shared
    r.tightinset(row,col,3)=r.tightinset(row,col,3)+hacbinset(3);
   end
  else
   error('Acessing orientation specification for local colorbar')
  end
 elseif ~isempty(hahcb) & ishandle(hahcb) & isempty(haorientation)
  error('Accessing colorbar orientation')
 end

 % determine space required if MATLAB legend or colorbar used with individual axes
 Tag=get(getappdata(ha,'LegendColorbarOuterList'),'Tag');
 if strcmpi(Tag,'Legend') | strcmpi(Tag,'Colorbar')
  if strcmpi(class(Tag),'cell')
   error('Only able to use a single MATLAB legend or colorbar per axes with ridgepack_multiplot')
  elseif ~strcmpi(get(getappdata(ha,'LegendColorbarOuterList'),'Units'),'normalized')
   error('Legend not in normalized coordinates')
  else
   r.legcolhandle(row,col)=getappdata(ha,'LegendColorbarOuterList');
  end
 end

 % position outside MATLAB legend or colorbar manually
 if isfield(r,'legcolhandle') && ishandle(r.legcolhandle(row,col))

  userdata=get(r.legcolhandle(row,col),'UserData');
  legpos=get(r.legcolhandle(row,col),'Position');
  legtight=get(r.legcolhandle(row,col),'TightInset');
  location=get(r.legcolhandle(row,col),'Location');

  if isfield(userdata,'Location') & any(strcmpi({'none','manual'},location))
   location=userdata.Location;
  elseif isempty(location)
   error('Unable to find legend location')
  else
   userdata.Location=location;
  end

  if ~isempty(strfind(location,'Outside'))
   if strcmpi(location,'NorthOutside')
    legpos(1)=position(1)+(position(3)-legpos(3))/2;
    legpos(2)=position(2)+position(4)+r.tightinset(row,col,4)+legtight(2);
    r.tightinset(row,col,4)=r.tightinset(row,col,4)+legtight(2)+legpos(4)+legtight(4);
    r.tightinset(row,col,1)=max(r.tightinset(row,col,1),legtight(1)+(legpos(3)-position(3))/2);
    r.tightinset(row,col,3)=max(r.tightinset(row,col,3),legtight(3)+(legpos(3)-position(3))/2);
   elseif strcmpi(location,'EastOutside') 
    legpos(1)=position(1)+position(3)+r.tightinset(row,col,3)+legtight(1);
    legpos(2)=position(2)+(position(4)-legpos(4))/2;
    r.tightinset(row,col,3)=r.tightinset(row,col,3)+legtight(1)+legpos(3)+legtight(3);
    r.tightinset(row,col,2)=max(r.tightinset(row,col,2),legtight(2)+(legpos(4)-position(4))/2);
    r.tightinset(row,col,4)=max(r.tightinset(row,col,4),legtight(4)+(legpos(4)-position(4))/2);
   elseif strcmpi(location,'WestOutside')
    legpos(1)=position(1)-r.tightinset(row,col,1)-legtight(3)-legpos(3);
    legpos(2)=position(2)+(position(4)-legpos(4))/2;
    r.tightinset(row,col,1)=r.tightinset(row,col,1)+legtight(1)+legpos(3)+legtight(3);
    r.tightinset(row,col,2)=max(r.tightinset(row,col,2),legtight(2)+(legpos(4)-position(4))/2);
    r.tightinset(row,col,4)=max(r.tightinset(row,col,4),legtight(4)+(legpos(4)-position(4))/2);
   elseif strcmpi(location,'SouthOutside') 
    legpos(1)=position(1)+(position(3)-legpos(3))/2;
    legpos(2)=position(2)-r.tightinset(row,col,2)-legtight(4)-legpos(4);
    r.tightinset(row,col,2)=r.tightinset(row,col,2)+legtight(2)+legpos(4)+legtight(4);
    r.tightinset(row,col,1)=max(r.tightinset(row,col,1),legtight(1)+(legpos(3)-position(3))/2);
    r.tightinset(row,col,3)=max(r.tightinset(row,col,3),legtight(3)+(legpos(3)-position(3))/2);
   else
    error('Only able to use North-, East-, South- or WestOutside with ridgepack_multiplot')
   end
   set(r.legcolhandle(row,col),'Position',legpos);
   set(r.legcolhandle(row,col),'UserData',userdata);
  end

 end

elseif ~isempty(ha) & isempty(r)

 error('Unable to access tightinset information')

end

if debug; disp(['...Leaving ',mfilename]); end


