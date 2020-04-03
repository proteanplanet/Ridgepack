
clf
clear

centlat=80; % degrees north
centlon=-150; % degrees east
horizon=60; % degrees of satellite horizon (0-90)

% location of grid file
gridloc='/Users/afroberts/data/MODEL/E3SM/DECK/grid';
gridfile='E3SM_LR_V1_grid.nc';
cgrid=true; 

% location of sea ice data
dataloc='/Users/afroberts/data/MODEL/E3SM/DECK/monthly/h1/archive/ice/reduced';
datafile='mpascice.hist.am.timeSeriesStatsMonthly.1980-03-01.nc';
varc='timeMonthly_avg_iceAreaCell';
datatitle='Sea Ice Extent';
threshold=0.15;

% plot location
plotloc='/Users/afroberts/work';

% %%%%%%%%%%%%%%% CHANGE ABOVE THIS LINE %%%%%%%%%%%%%%%%%% %

% obtain grid information (vertices, cell centers)
cd(gridloc)
ncvert=ridgepack_clone(gridfile,{'latVertex','lonVertex','dcEdge',...
                                 'verticesOnCell','indexToCellID',...
                                 'nEdgesOnCell','edgesOnCell',...
                                 'cellsOnEdge'});

nccell=ridgepack_clone(gridfile,{'latCell','lonCell'});

% obtain data
cd(dataloc)
ncc=ridgepack_clone(datafile,varc);

% put up satview
ridgepack_satview(centlat,centlon,horizon,1,2);

% plot ice edge 
ridgepack_pthresholde3sm(ncc,varc,threshold,ncvert,...
                         centlat,centlon,horizon)

% plot coastal outline
ridgepack_psatcoaste3sm(ncvert,cgrid,...
                        centlat,centlon,horizon);

% add title
title([datatitle])

cd(plotloc)
ridgepack_fprint('png','Outfile_Sea_Ice_Threshold',1,1)


