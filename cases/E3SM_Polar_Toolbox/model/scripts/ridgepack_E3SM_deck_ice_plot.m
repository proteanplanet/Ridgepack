
clf
clear

% Arctic Center
centlat=80; % degrees north
centlon=-150; % degrees east
horizon=30; % degrees of satellite horizon (0-90)
altitude=1; % Mean Earth radius multiple
cgrid=true; % plot c-grid coastline
coastname='DECK'; % grid name

% location of grid file
gridloc='/Users/afroberts/data/MODEL/E3SM/DECK/grid';
gridfile='E3SM_LR_V1_grid.nc';

% location of sea ice data
dataloc='/Users/afroberts/data/MODEL/E3SM/DECK/monthly/h1/archive/ice/reduced';
datafile='mpascice.hist.am.timeSeriesStatsMonthly.1980-03-01.nc';
choice=3;

% choices for sea ice model
if choice==1
 field='timeMonthly_avg_iceAreaCell';
 cont=[0:0.05:1];
 ref=0.5;
 units='';
 datatitle='Sea Ice Area';
elseif choice==2
 field='timeMonthly_avg_iceVolumeCell';
 cont=[0:0.25:5.25];
 ref=0;
 units='m';
 datatitle='Mean Grid-Cell Sea Ice Thickness';
elseif choice==3
 field='timeMonthly_avg_divergence';
 cont=[-1.0:0.1:-0.1 0.1:0.1:1];
 ref=0;
 units='day^{-1}';
 datatitle='Sea Ice Divergence';
else
 error('data choice not available')
end

% plot location
plotloc='/Users/afroberts/work';

% %%%%%%%%%%%%%%% CHANGE ABOVE THIS LINE %%%%%%%%%%%%%%%%%% %

% obtain grid information
cd(gridloc)
ncvert=ridgepack_clone(gridfile,{'latVertex','lonVertex','dcEdge',...
                                 'verticesOnCell','indexToCellID',...
                                 'nEdgesOnCell','edgesOnCell',...
                                 'cellsOnEdge'});

% obtain field data
cd(dataloc)
ncdata=ridgepack_clone(datafile,field);

% set up satellite view
ridgepack_satview(centlat,centlon,horizon,1,2);

% plot cell resolution
ridgepack_psatcole3sm(ncdata,field,ncvert,cont,ref,...
                      centlat,centlon,horizon,altitude,false);

% plot coastal outline
ridgepack_psatcoaste3sm(ncvert,cgrid,coastname,...
                             centlat,centlon,horizon);

% add colorbar
ridgepack_colorbar(cont,units);

title([datatitle,' (',ncdata.(field).long_name,')'])

cd(plotloc)
ridgepack_fprint('png','Outfile',1,1)


