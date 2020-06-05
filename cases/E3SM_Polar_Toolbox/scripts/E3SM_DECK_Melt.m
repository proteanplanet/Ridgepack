
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

ridgepack_polarm('antarctic','noland','grid','label');
geoshow(SCP,'FaceColor',0.95*[1 1 1],'EdgeColor',0.5*[1 1 1])
mask=find(nc.timeMonthly_avg_iceAreaCell.data>0.001);
ridgepack_e3smcolors(nc,...
                  'timeMonthly_avg_freezingMeltingPotential',...
                   ncvert,mask,[-90 -60 -30:2:4],'linear',-14);
h=geoshow(STC,'Color','m','LineWidth',1)
legend(h,'Extent','Location','SouthEast')
legend boxoff
title('September DECK PI Years 0400-0499 Mean Freeze-Melt Potential')
ridgepack_fprint('png','Antarctic_Sept_0400_0499_PI_DECK_FMP',1,2)

clf

ridgepack_polarm('antarctic','noland','grid','label');
geoshow(SCP,'FaceColor',0.95*[1 1 1],'EdgeColor',0.5*[1 1 1])
mask=find(nc.timeMonthly_avg_iceAreaCell.data>0.001);
ridgepack_e3smcolors(nc,...
                  'timeMonthly_avg_basalIceMelt',...
                   ncvert,mask,[0 0.01 0.03 0.06 0.12 0.18:0.18:1.80+0.18]*10^-7);
h=geoshow(STC,'Color','m','LineWidth',1)
legend(h,'Extent','Location','SouthEast')
legend boxoff
title('September DECK PI Years 0400-0499 Mean Basal Melt')
ridgepack_fprint('png','Antarctic_Sept_0400_0499_PI_DECK_Basal_Melt',1,2)

clf


ridgepack_polarm('antarctic','noland','grid','label');
geoshow(SCP,'FaceColor',0.95*[1 1 1],'EdgeColor',0.5*[1 1 1])
mask=find(nc.timeMonthly_avg_iceAreaCell.data>0.001);
ridgepack_e3smcolors(nc,...
                  'timeMonthly_avg_lateralIceMelt',...
                   ncvert,mask,[0 0.001 0.01 0.025 0.05 0.075 0.1:0.1:1.1]*10^-9);
h=geoshow(STC,'Color','m','LineWidth',1)
legend(h,'Extent','Location','SouthEast')
legend boxoff
title('September DECK PI Years 0400-0499 Mean Lateral Melt')
ridgepack_fprint('png','Antarctic_Sept_0400_0499_PI_DECK_Lateral_Melt',1,2)

clf

ridgepack_polarm('antarctic','noland','grid','label');
geoshow(SCP,'FaceColor',0.95*[1 1 1],'EdgeColor',0.5*[1 1 1])
mask=find(nc.timeMonthly_avg_iceAreaCell.data>0.001);
ridgepack_e3smcolors(nc,...
                  'timeMonthly_avg_snowMelt',...
                   ncvert,mask,[0 0.5 1:1:13]*10^-9);
h=geoshow(STC,'Color','m','LineWidth',1)
legend(h,'Extent','Location','SouthEast')
legend boxoff
title('September DECK PI Years 0400-0499 Mean Snow Melt')
ridgepack_fprint('png','Antarctic_Sept_0400_0499_PI_DECK_Snow_Melt',1,2)

clf

ridgepack_polarm('antarctic','noland','grid','label');
geoshow(SCP,'FaceColor',0.95*[1 1 1],'EdgeColor',0.5*[1 1 1])
mask=find(nc.timeMonthly_avg_iceAreaCell.data>0.001);
ridgepack_e3smcolors(nc,...
                  'timeMonthly_avg_surfaceIceMelt',...
                   ncvert,mask,[0 0.001 0.01 0.1:0.1:1.1]*10^-9);
h=geoshow(STC,'Color','m','LineWidth',1)
legend(h,'Extent','Location','SouthEast')
legend boxoff
title('September DECK PI Years 0400-0499 Mean Surface Melt')
ridgepack_fprint('png','Antarctic_Sept_0400_0499_PI_DECK_Surface_Melt',1,2)

clf

ridgepack_polarm('antarctic','noland','grid','label');
geoshow(SCP,'FaceColor',0.95*[1 1 1],'EdgeColor',0.5*[1 1 1])
mask=find(nc.timeMonthly_avg_iceAreaCell.data>0.001);
ridgepack_e3smcolors(nc,...
                  'timeMonthly_avg_frazilFormation',...
                   ncvert,mask,[0:0.1:1.5]*10^-8);
h=geoshow(STC,'Color','m','LineWidth',1)
legend(h,'Extent','Location','SouthEast')
legend boxoff
title('September DECK PI Years 0400-0499 Mean Frazil Formation')
ridgepack_fprint('png','Antarctic_Sept_0400_0499_PI_DECK_Frazil_Formation',1,2)

clf

ridgepack_polarm('antarctic','noland','grid','label');
geoshow(SCP,'FaceColor',0.95*[1 1 1],'EdgeColor',0.5*[1 1 1])
mask=find(nc.timeMonthly_avg_iceAreaCell.data>=0.0);
ridgepack_e3smcolors(nc,...
                  'timeMonthly_avg_seaSurfaceSalinity',...
                   ncvert,mask,[32.3:0.1:34.4]);
h=geoshow(STC,'Color','m','LineWidth',1)
legend(h,'Extent','Location','SouthEast')
legend boxoff
title('September DECK PI Years 0400-0499 Mean SSS')
ridgepack_fprint('png','Antarctic_Sept_0400_0499_PI_DECK_SSS',1,2)

clf

ridgepack_polarm('antarctic','noland','grid','label');
geoshow(SCP,'FaceColor',0.95*[1 1 1],'EdgeColor',0.5*[1 1 1])
mask=find(nc.timeMonthly_avg_iceAreaCell.data>=0.0);
ridgepack_e3smcolors(nc,...
                  'timeMonthly_avg_seaSurfaceTemperature',...
                   ncvert,mask,[-1.8 -1.5:0.5:0 1:1:11]);
h=geoshow(STC,'Color','m','LineWidth',1)
legend(h,'Extent','Location','SouthEast')
legend boxoff
title('September DECK PI Years 0400-0499 Mean SST')
ridgepack_fprint('png','Antarctic_Sept_0400_0499_PI_DECK_SST',1,2)

clf

ridgepack_polarm('antarctic','noland','grid','label');
geoshow(SCP,'FaceColor',0.95*[1 1 1],'EdgeColor',0.5*[1 1 1])
mask=find(nc.timeMonthly_avg_iceAreaCell.data>0.001);
ridgepack_e3smcolors(nc,...
                  'timeMonthly_avg_surfaceTemperatureCell',...
                   ncvert,mask,[-30:2:0]);
h=geoshow(STC,'Color','m','LineWidth',1)
legend(h,'Extent','Location','SouthEast')
legend boxoff
title('September DECK PI Years 0400-0499 Mean Surface Temp')
ridgepack_fprint('png','Antarctic_Sept_0400_0499_PI_DECK_SurfaceT',1,2)

clf

ridgepack_polarm('antarctic','noland','grid','label');
geoshow(SCP,'FaceColor',0.95*[1 1 1],'EdgeColor',0.5*[1 1 1])
mask=find(nc.timeMonthly_avg_iceAreaCell.data>0.001);
ridgepack_e3smcolors(nc,...
                  'timeMonthly_avg_oceanHeatFlux',...
                   ncvert,mask,[-140:10:-10 -5 -2 -1 0],'linear',0);
h=geoshow(STC,'Color','m','LineWidth',1)
legend(h,'Extent','Location','SouthEast')
legend boxoff
title('September DECK PI Years 0400-0499 Mean Ocean Heat Flux')
ridgepack_fprint('png','Antarctic_Sept_0400_0499_PI_DECK_OHF',1,2)

clf

ridgepack_polarm('antarctic','noland','grid','label');
geoshow(SCP,'FaceColor',0.95*[1 1 1],'EdgeColor',0.5*[1 1 1])
mask=find(nc.timeMonthly_avg_iceAreaCell.data>0.001);
ridgepack_e3smcolors(nc,...
                  'timeMonthly_avg_congelation',...
                   ncvert,mask,[0:0.1:1.1]*10^-7);
h=geoshow(STC,'Color','m','LineWidth',1)
legend(h,'Extent','Location','SouthEast')
legend boxoff
title('September DECK PI Years 0400-0499 Mean Congelation Growth')
ridgepack_fprint('png','Antarctic_Sept_0400_0499_PI_DECK_Congelation',1,2)

clf

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

clf

ridgepack_polarm('antarctic','noland','grid','label');
geoshow(SCP,'FaceColor',0.95*[1 1 1],'EdgeColor',0.5*[1 1 1])
mask=find(nc.timeMonthly_avg_iceAreaCell.data>0.15);
ridgepack_e3smcolors(nc,...
                  'timeMonthly_avg_iceVolumeCell',...
                   ncvert,mask,[0:0.25:3.25]);
h=geoshow(STC,'Color','m','LineWidth',1)
legend(h,'Extent','Location','SouthEast')
legend boxoff
title('September DECK PI Years 0400-0499 Mean Thickness')
ridgepack_fprint('png','Antarctic_Sept_0400_0499_PI_DECK_Thickness',1,2)





