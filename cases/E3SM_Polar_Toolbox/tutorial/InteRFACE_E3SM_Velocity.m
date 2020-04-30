% InteRFACE_E3SM_velocity: Tutorial Script
% 
% This script demonstrates how to quickly look at sea ice drift
% from E3SM, and output streamlines in a KML file.
%
% InteRFACE Tutorial, April 2020
% Andrew Roberts, LANL, afroberts@lanl.gov

clf
clear 

% CHANGE THESE THINGS TO MAKE THIS WORK FOR YOU: 

% plot location
plotloc='/Users/afroberts/work/tutorial';

% grab out the grid data
gridloc='/Users/afroberts/data/MODEL/E3SM/DECK/grid';
gridfile='E3SM_LR_V1_grid.nc';


% data file location
dataloc='/Users/afroberts/data/MODEL/E3SM/DECK/monthly/h1/archive/ice/hist';
datafile='mpascice.hist.am.timeSeriesStatsMonthly.2000-03-01.nc';

% EVERYTHING BELOW THIS LINE SHOULD JUST WORK

% get the mesh data
cd(gridloc)
ncvert=ridgepack_clone(gridfile,{'latVertex','lonVertex',...
                                'verticesOnCell','indexToCellID',...
                                'nEdgesOnCell','edgesOnCell',...
                                'cellsOnEdge','cellsOnVertex',...
                                'edgesOnVertex'});

nccell=ridgepack_clone(gridfile,{'latCell','lonCell'});

ncedge=ridgepack_clone(gridfile,{'latEdge','lonEdge'});

ncvert.latCell=nccell.latitude;
ncvert.lonCell=nccell.longitude;
ncvert.latEdge=ncedge.latitude;
ncvert.lonEdge=ncedge.longitude;

% grab the data
cd(dataloc)

% fields we would like
fields={'timeMonthly_avg_iceAreaCell',...
        'timeMonthly_avg_uVelocityGeo',...
        'timeMonthly_avg_vVelocityGeo'}
 
% ingest the data, and remove the time dimension immediately
nc=ridgepack_reduce(ridgepack_clone(datafile,fields),{'time'});

mask=find(ncvert.latCell.data<pi/2);

% grab the coastal data and 15% sea ice concentration contour
[nccoast,SCP,SCL,STP,STC]=ridgepack_e3smseasaw(ncvert,[nc],...
                       'timeMonthly_avg_iceAreaCell',...
                        0.15);

% now add streamlines
hemisphere=1;
STS=ridgepack_e3smstream(nc,'timeMonthly_avg_uVelocityGeo',...
                     nc,'timeMonthly_avg_vVelocityGeo',...
                     nc,'timeMonthly_avg_iceAreaCell',...
                     ncvert,hemisphere,7,[0:2:32]*10^-2,0)

% add in the 80% contour
h=geoshow(STC,'Color',[0.8 0.33 0],'LineWidth',1)

% add in the native model coastline
geoshow(SCL,'Color',0.5*[1 1 1],'LineWidth',1)

% add a legend
legend(h,'15% concentration','Location','SouthEast')
legend boxoff

% add the title
title('March 2000 Sea Ice Drift, E3SM DECK h1')

% plot data
cd(plotloc)
ridgepack_fprint('png',mfilename,1,2)

% write a KMLfile of streamlines
kmlwrite('Drift',STS,'Color','c','LineWidth',1, ...
         'Description','Sea Ice Drift','Name','MPAS Sea Ice');

% write a netcdf coastal file
ridgepack_write(nccoast,'E3SM_Coast')




