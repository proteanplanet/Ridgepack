function ridgepack_title(nc,text,pos)

% ridgepack_title - Add a title to a plot made in ncbox
%
% function ridgepack_title(nc,text,pos)
%
% Input:
% nc   - netcdf structure (see ridgepack_struct for more information)
% text - text to be added to the title
% pos  - position of title (1,2 or 3)
%        1 - top rung
%        2 - bottom rung
%        3 - remove previous title
%
% Output:
% Output is graphical only
%
% Andrew Roberts, Naval Postgraduate School, March 2018  (afrobert@nps.edu)
%

global debug;
if debug; disp(['Entering ',mfilename,'...']); end

if ishold
	currenttitle=get(get(gca,'Title'),'String');
else
	currenttitle='';
end

% clear title for identical fields overplotted
repeatfields={'Shading:','Contour:','Vectors:'};
for i=1:length(currenttitle)

 if ~iscell(currenttitle) ; 
   error('Unable to access current title - turn hold off'); 
 end

 name=char(currenttitle{i});
 if ~isempty(name) & length(name)>=8
  if any(strcmp(text(1:8),repeatfields)) & strcmp(name(1:8),text(1:8))
   currenttitle{i}='';
  end
 end

end

if pos==1
	notpos=2;
elseif pos==2
	notpos=1;
elseif pos==3
	pos=1;
	notpos=2;
	currenttitle='';
end

if iscell(currenttitle)

 if length(currenttitle)==2

	if isempty(char(currenttitle{pos}))
		 currenttitle{pos}=text;
	elseif not(strcmp(text,char(currenttitle{pos})))
		 currenttitle{notpos}=currenttitle{pos};
		 currenttitle{pos}=text;
	end

	h=title(currenttitle,'HorizontalAlignment','center','FontWeight','normal');

 elseif length(currenttitle)==1

	newtitle{pos}=text;
	newtitle{notpos}=char(currenttitle);

	h=title(newtitle,'HorizontalAlignment','center','FontWeight','normal');

 else

	newtitle{pos}=text;
	newtitle{notpos}='';

	h=title(newtitle,'HorizontalAlignment','center','FontWeight','normal');

 end

else

	newtitle{pos}=text;
	newtitle{notpos}='';

	h=title(newtitle,'HorizontalAlignment','center','FontWeight','normal');

end

if debug; disp(['...Leaving ',mfilename]); end

