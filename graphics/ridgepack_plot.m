function ridgepack_plot(nc,X,Y,dimnames,bounds,options)

% ridgepack_plot - Graphs data from an nc structure on Cartesian axes
%
% function ridgepack_plot(nc,X,Y,dimnames,bounds,options)
%
% INPUT:
%
% nc       - netcdf data structure
%
% X        - character variable of x variable in nc structure
%
% Y        - character variable of y variable in nc structure
%
% dimnames - cell array of the dimensions to be removed from the nc
%            structure to obtain the requested cross section. See example
%            provided under ridgepack_reduce for more information.
%
% bounds   - cell array of the boundaries on each dimension reduction.
%            This includes either a low and high index over which a mean
%            or slice is to be taken, or a low and high value representing
%            the range over which the data is to be extracted.  Time
%            values can be used. See example provided for ridgepack_reduce for
%            more information.
%
% options  - plot options used in Matlab function "plot".
%            This variable is option.
% 
%
% OUTPUT:
%
% Output is graphical
%
% Ridgepack Version 1.0
% Andrew Roberts, Naval Postgraduate School, March 2018 (afrobert@nps.edu)
% 

global debug;
if debug; disp(['Entering ',mfilename,'...']); end

if ~isstruct(nc)
	error('Not a netcdf structure')
end

if nargin<4; dimnames={}; end
if nargin<5; bounds={}; end
if nargin<6; options=''; end

% generate data
[nc]=ridgepack_select(nc,X,Y,Y,dimnames,bounds);

% get all previous legend entries if a hold is in place
i=0;
legendstring=cell(1,1);
if ishold 
	disp('Getting previous legend information')
	hhh=get(gca);
	legendstring=cell(length(hhh.Children),1);
	for i=1:length(hhh.Children)
		hgg=get(hhh.Children(i));
		legendstring{i}=hgg.DisplayName;
	end
end

if any(size(nc.(Y).data)==1)
 legendstring{i+1}=nc.attributes.title;
else
 legendstring{i+1}=nc.(Y).long_name;
end

% change units if need be
[x,nc]=ridgepack_standardunits(nc,X);
[y,nc]=ridgepack_standardunits(nc,Y);

% plot the data
plot(nc.(X).data,nc.(Y).data,options);

% add legend, whiting out the legend color if only one trace is shown
disp('Adding legend')
if ishold
	legend(legendstring,'Location','NorthOutside');
	legend('boxoff');
elseif size(nc.(Y).data,1)>1 && size(nc.(Y).data,2)>1 
	title(legendstring);
	if size(nc.(Y).data,1)==length(nc.(X).data)
		units=ridgepack_units(nc,char(nc.(Y).dimension{2}));
		for i=1:size(nc.(Y).data,2)
			leg{i}=[num2str(nc.(char(nc.(Y).dimension{2})).data(i))];
		end
		leg{1}=[char(leg(1)),'  ',units];
		legend(leg,'Location','northeast');
	elseif size(nc.(Y).data,2)==length(nc.(X).data)
		units=ridgepack_units(nc,char(nc.(Y).dimension{1}));
		for i=1:size(nc.(Y).data,1)
			leg{i}=[num2str(nc.(char(nc.(Y).dimension{1})).data(i))];
		end
		leg{1}=[char(leg(1)),'  ',units];
		legend(leg,'Location','northeast');
	else
		disp('Unable to work out the legend');
	end
else
	[legend_h,object_h]=legend(legendstring,'Location','NorthOutside');
	legend('boxoff');
end

% label axes
ridgepack_labelxy(nc,X,Y,true);

drawnow

if debug; disp(['...Leaving ',mfilename]); end


