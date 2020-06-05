
clf
clear 

cd('/Users/afroberts/data/MODEL/E3SM/DECK/monthly/PI/archive/ice/hist') 

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


[nccoast,SCP,STP,STC]=ridgepack_e3smseasaw(ncvert,[nc],...
                       'timeMonthly_avg_iceAreaCell',...
                        0.15);

plotloc='/Users/afroberts/work';
cd(plotloc)

if 1==0

ridgepack_polarm('antarctic','noland','grid','label');
geoshow(SCP,'FaceColor',0.95*[1 1 1],'EdgeColor',0.5*[1 1 1])
mask=find(nc.timeMonthly_avg_iceAreaCell.data>0.15);
ridgepack_e3smcolors(nc,...
                  'timeMonthly_avg_iceAreaCell',...
                   ncvert,mask,[0:0.05:1]);
h=geoshow(STC,'Color','m','LineWidth',1)
legend(h,'Extent','Location','SouthEast')
legend boxoff
title('September DECK PI Years 0400-0499 Mean Concentration')
ridgepack_fprint('png','Antarctic_Sept_0400_0499_PI_DECK_Concentration',1,2)

end

clf

ridgepack_polarm('antarctic','noland','grid','label');
geoshow(SCP,'FaceColor',0.95*[1 1 1],'EdgeColor',0.5*[1 1 1])
mask=find(nc.timeMonthly_avg_iceAreaCell.data>0.15);
ridgepack_e3smcolors(nc,...
                  'timeMonthly_avg_iceVolumeCell',...
                   ncvert,mask,[0:0.25:3.25]);
mask=find(nc.timeMonthly_avg_iceAreaCell.data<0.85);
ridgepack_e3smeshs(ncvert,mask,'w')
h=geoshow(STC,'Color','m','LineWidth',1);
legend(h,'Extent','Location','SouthEast')
legend boxoff
title('September DECK PI Years 0400-0499 Mean Thickness')
ridgepack_fprint('png','Antarctic_Sept_0400_0499_PI_DECK_Thickness',1,2)





