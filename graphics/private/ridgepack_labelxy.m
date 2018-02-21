function ridgepack_labelxy(nc,X,Y,graph)

% function ridgepack_labelxy(nc,X,Y,graph)
%
% This function is part of Ridgepack Version 1.0.
% It labels figure axes
% 
% INPUT:
%
% nc    - netcdf structure (see ncstruct for more details)
% X     - x-coordinate of Cartesian plot in nc
% Y     - y-coordinate of Cartesian plot in nc
% graph - logical which is true if a graph is being plotted
%         or else may be omitted.
%
% OUTPUT:
%
% All output is graphical and appears on the current selected axes
%
% Andrew Roberts, Naval Postgraduate School, March 2018 (afrobert@nps.edu)

% set defaults
if nargin<4
	graph=false;
end

% derive the best font size
figout=get(gca,'OuterPosition');
fontsize=min(10,max(6,10*sqrt((figout(3).^2)+(figout(4).^2))/sqrt(2)));

% label axes and check set boundaries
if (isfield(nc.(X),'units') && ~isempty(nc.(X).units)) | strcmp(X,'time')
 units=ncunits(nc,X);
 if strcmp(X,'time')
       	xlabel('Time','FontSize',fontsize)
	set(gca,'xlim',[min(nc.(X).data(:)) max(nc.(X).data(:))])
	%datetick('x','keeplimits')
	datetick('x','keepticks','keeplimits')
	pos=get(gca,'Position');
	out=get(gca,'OuterPosition');

	if ~graph
	 % cover up power-of-ten symbol on plot
	 patch([pos(1) ; pos(1) ; pos(1)+pos(3) ; pos(1)+pos(3)],...
	       [pos(2) ; out(2) ; out(2) ; pos(2)],...
	       [1 1 1],'linestyle','none') ;
        end

 elseif strcmp(char(nc.(X).units),'degrees_north')
       	xlabel('Latitude (^{o}N)','FontSize',fontsize)
 elseif strcmp(char(nc.(X).units),'degrees_east')
       	xlabel('Longitude (^{o}E)','FontSize',fontsize)
 elseif ~strcmp(units,'')
       	xlabel([char(nc.(X).long_name),' (',units,')'],'FontSize',fontsize)
 else
 	xlabel([char(nc.(X).long_name)],'FontSize',fontsize)
 end
else
 xlabel([char(nc.(X).long_name)],'FontSize',fontsize)
end

if graph & ~strcmp(Y,'time') & ...
   ~isempty(findstr(Y,'latitude')) & ~isempty(findstr(Y,'longitude'))
	if isfield(nc.(Y),'units') && ~isempty(nc.(Y).units)
		units=ncunits(nc,Y);
		ylabel(units,'FontSize',fontsize)
	end
elseif (isfield(nc.(Y),'units') && ~isempty(nc.(Y).units)) | strcmp(Y,'time')
 units=ncunits(nc,Y);
 if strcmp(Y,'time')
       	ylabel('Time','FontSize',fontsize)
	set(gca,'ylim',[min(nc.(Y).data(:)),max(nc.(Y).data(:))])
	datetick('y','keeplimits')
	pos=get(gca,'Position');
	out=get(gca,'OuterPosition');

	if ~graph
	 % cover up power-of-ten symbol on plot
	 patch([pos(1) ; pos(1) ; pos(1)+pos(3) ; pos(1)+pos(3)],...
	       [pos(2)+pos(4) ; out(2)+out(4) ; out(2)+out(4) ; pos(2)+pos(4)],...
	       [1 1 1],'linestyle','none') ;
        end

 elseif strcmp(char(nc.(Y).units),'degrees_north')
       	ylabel('Latitude (^{o}N)','FontSize',fontsize)
 elseif strcmp(char(nc.(Y).units),'degrees_east')
       	ylabel('Longitude (^{o}E)','FontSize',fontsize)
 elseif ~strcmp(units,'')
	ylabel([char(nc.(Y).long_name),' (',units,')'],'FontSize',fontsize)
 else
 	ylabel([char(nc.(Y).long_name)],'FontSize',fontsize)
 end
else
 ylabel([char(nc.(Y).long_name)],'FontSize',fontsize)
end

