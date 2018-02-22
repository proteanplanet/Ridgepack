function ridgepack_multiplot(nrows,ncols,row,col,notation)

% ridgepack_multiplot - Efficiently plot multiple axes, colorbars and legends on one figure
%
% function ridgepack_multiplot(nrows,ncols,row,col,notation)
%
% This replaces much of the functionality of the MATLAB function 'subplot'. It plots
% aligned axes in rows and columns on the current figure, and adds a notation inside 
% the top left corner of each axes. Axes are progressively aligned and packed as they 
% are added to the figure, with space allocated for colorbars, axes labeling, and 
% titles.  When map axes are plotted, they must have the same width and height 
% proportions as their top (column) or right (row) axes, respectively. Axes of 
% different proportions may be plotted along the top row and left column. When
% plotting of a figure is complete, the function ridgepack_multialign completes alignment. 
% The function ridgepack_multialign can also add a title to all plots.  A colorbar 
% associated with a particular axes can be switched into a colorbar that applies
% to the entire multiplot. To do this, use the function ridgepack_multicb. 
%
% Inputs:
% nrows    - total number of rows of axes in the figure 
% ncols    - total number of columns of axes in the figure
% row      - row number being added to the figure
% col      - column number being added to the figure
% notation - string to be added to the top left of the axes as notation (optional)
%
% Example: 
% ridgepack_multiplot(1,2,1,1,'a') , ridgepack_polarm('seaice');
% ridgepack_multiplot(1,2,1,2,'b') , ridgepack_polarm('seaice');
% ridgepack_multialign(gcf,'This is the title')
% This will plot a figure with two aligned maps, an 'a' in the top left 
% of the left-hand axes, and a 'b' in the top left of the righ-hand axes.
% Once all calls have been made to ridgepack_multiplot, a final call of ridgepack_multialign 
% brings them into alignement and uses all available space on the figure. 
%
% When to use ridgepack_multiplot and when to use subplot:
% The difference between ridgepack_multiplot and subplot is that the latter does not
% conduct any fitting of axes on the figure, whereas ridgepack_multiplot fits axes
% together after they have been plotted, removing the laborious work of drawing
% publication-ready figures or model output with many frames per figure.  
% Subplot functionality does allow for axes to span two or more rows and columns,
% whereas ridgepack_multiplot does not, and this is where subplot should be used.
% Finally, ridgepack_multiplot has been set up to explicitly handle the icepack colorbar,
% ensuring correct placement and spacing around colorbars, regardless of where they 
% are positioned in the figure window.
%
% Special case of single column multiple axes figures:
% If ncols=1 (a single column figure), then the default span of the axes
% is to stretch across the entire figure, as is useful for comparing graphs
% stacked upon each other.
%
% Order of adding axes to a figure:
% Each time this function is run with row=1 and col=1, the underlying 
% multiplot configuration is reset.  As a result, you must always call
% this function for row=1 and col=1 first, before adding other axes
% to the plot.  All other axes can be added in any order.
%
% Ridgepack Version 1.0
% Andrew Roberts, Naval Postgraduate School, March 2018 (afrobert@nps.edu)
%

global debug;
if debug; disp(['Entering ',mfilename,'...']); end

if nargin<4
 error('Must provide four arguments: nrows, ncols, row, col')
elseif ~isnumeric(nrows) | nrows<1
 error('nrows must be a positive integer')
elseif ~isnumeric(ncols) | ncols<1
 error('ncolss must be a positive integer')
elseif ~isnumeric(row) | row<1 | row>nrows
 error(['row must be a positive integer less than or equal to ',num2str(nrows)])
elseif ~isnumeric(col) | col<1 | col>ncols
 error(['col must be a positive integer less than or equal to ',num2str(ncols)])
else
 maxrc=max(nrows,ncols);
end

if nargin<5
 notation='';
elseif ~ischar(notation)
 error('Notation should be entered as a string')
end

if nargin<6
 xtendcol=0;
elseif ~isnumeric(xtendcol) | xtendcol>1 | xtendcol<0
 error('xtendcol must be either 0 or 1')
elseif col>1
 error('Can only extend axes for col=1')
end

if row==1 & col==1
 r=getappdata(gcf,'MultiplotConfiguration');
 if ~isempty(r); rmappdata(gcf,'MultiplotConfiguration'); end
end

ha=ridgepack_multipos(gcf,[],nrows,ncols,row,col,notation);
hold(ha,'on') 

hlf=getappdata(gcf,'MultiplotFigureListener');
if isempty(hlf)
 hlf=addlistener(gcf,{'CurrentAxes'},'PostSet',@ridgepack_multifix);
 setappdata(gcf,'MultiplotFigureListener',hlf)
end

if debug; 
 disp(['...Leaving ',mfilename]); 
end

