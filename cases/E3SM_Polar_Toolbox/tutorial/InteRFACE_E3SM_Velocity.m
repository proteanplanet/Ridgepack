% InteRFACE_E3SM_Velocity: Tutorial Script
% 
% This script demonstrates how to quickly look at sea ice drift
% from E3SM.
%
% InteRFACE Tutorial, April 2020
% Andrew Roberts, LANL, afroberts@lanl.gov

clf
clear 



% grab out the grid data
gridloc='/Users/afroberts/data/MODEL/E3SM/DECK/grid';
gridfile='E3SM_LR_V1_grid.nc';




cd(gridloc)
ncvert=ridgepack_clone(gridfile,{'latVertex','lonVertex',...
                                'verticesOnCell','indexToCellID',...
                                'nEdgesOnCell','edgesOnCell',...
                                'cellsOnEdge','cellsOnVertex',...
                                'edgesOnVertex'});

nccell=ridgepack_clone(gridfile,{'latCell','lonCell',...
                                'verticesOnCell','indexToCellID',...
                                'nEdgesOnCell','edgesOnCell',...
                                'cellsOnEdge','cellsOnVertex',...
                                'edgesOnVertex'});

ncedge=ridgepack_clone(gridfile,{'latEdge','lonEdge'});

ncvert.latCell=nccell.latitude;
ncvert.lonCell=nccell.longitude;
ncvert.latEdge=ncedge.latitude;
ncvert.lonEdge=ncedge.longitude;

nccell.latVertex=ncvert.latitude;
nccell.lonVertex=ncvert.longitude;
nccell.latEdge=ncedge.latitude;
nccell.lonEdge=ncedge.longitude;

% grab data
dataloc='/Users/afroberts/data/MODEL/E3SM/DECK/monthly/PI/archive/ice/hist';
datafile='mpascice.hist.am.DECK.melting.mean.09.0400_0499.nc';
cd(dataloc)
nc=ridgepack_clone(datafile);

mask=find(nccell.latitude.data<pi/2);

[nccoast,SCP,STP,STC]=ridgepack_e3smseasaw(ncvert,[nc],...
                       'timeMonthly_avg_iceAreaCell',...
                        0.15);

ridgepack_polarm('antarctic','noland','grid','label');

geoshow(SCP,'FaceColor',0.95*[1 1 1],'EdgeColor',0.5*[1 1 1])

ridgepack_e3smstream(nc,'timeMonthly_avg_uVelocityGeo',...
                     nc,'timeMonthly_avg_vVelocityGeo',...
                     nc,'timeMonthly_avg_iceAreaCell',...
                     ncvert,hemisphere,7,[0:2:32]*10^-2,0)

h=geoshow(STC,'Color',[0.8 0.33 0],'LineWidth',1)

legend(h,'Extent','Location','SouthEast')
legend boxoff

title('September DECK PI Years 0400-0499 Mean Sea Ice Drift')

% plot data
plotloc='/Users/afroberts/work';
cd(plotloc)
ridgepack_fprint('png','Antarctic_Sept_0400_0499_PI_DECK_Drift',1,2)





