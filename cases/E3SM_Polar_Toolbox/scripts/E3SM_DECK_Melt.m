
clf
clear 

% grab out the grid data
gridloc='/Users/afroberts/data/MODEL/E3SM/highres/grid';
gridfile='mpaso.rst.nc';

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
dataloc='/Users/afroberts/Work';

timeslice=31;

timetype='Daily';

titletag='E3SM-HR';

projtag='antarctic';

datafile=['mpascice.hist.am.timeSeriesStats',...
          timetype,'.1950-01-01.nc'];

cd(dataloc)

iceareacell=['time',timetype,'_avg_iceAreaCell'];

nca=ridgepack_reduce(ridgepack_clone(datafile,...
                     {iceareacell},timeslice),{'time'});

nct=ridgepack_clone(datafile,...
                    {'timeDaily_counter',...
                     'xtime_startDaily',...
                     'xtime_endDaily'});

stime=nct.xtime_startDaily.data(timeslice,1:10);
timestamp=[stime,'_',timetype];

[nccoast,SCP,SCL,STP,STC]=ridgepack_e3smseasaw(ncvert,[nca],...
                       iceareacell,0.15);

plotloc='/Users/afroberts/work';
cd(plotloc)

ridgepack_polarm(projtag,'noland','grid','label');
geoshow(SCP,'FaceColor',0.95*[1 1 1],'EdgeColor',0.5*[1 1 1])
mask=find(nca.(iceareacell).data>0.001);
field=['time',timetype,'_avg_freezingMeltingPotential'];
nc=ridgepack_reduce(ridgepack_clone(datafile,...
                    {field},timeslice),{'time'});
ridgepack_e3smcolors(nc,field,ncvert,mask,...
                     [-90 -60 -30:2:4],'linear',-14);
h=geoshow(STC,'Color','m','LineWidth',1);
legend(h,'Extent','Location','SouthEast')
legend boxoff
titlefield='Freeze-Melt Potential';
filetag='FreezeMeltPotential';
title([titletag,' ',stime,' ',timetype,' ',titlefield])
ridgepack_fprint('png',[projtag,'_',titletag,'_',...
                        stime,'_',timetype,'_',filetag],1,2)

clf

ridgepack_polarm(projtag,'noland','grid','label');
geoshow(SCP,'FaceColor',0.95*[1 1 1],'EdgeColor',0.5*[1 1 1])
mask=find(nca.(iceareacell).data>0.001);
field=['time',timetype,'_avg_basalIceMelt'];
nc=ridgepack_reduce(ridgepack_clone(datafile,...
                    {field},timeslice),{'time'});
ridgepack_e3smcolors(nc,field,ncvert,mask,...
                   [0 0.01 0.03 0.06 0.12 0.18:0.18:1.80+0.18]*10^-7);
h=geoshow(STC,'Color','m','LineWidth',1)
legend(h,'Extent','Location','SouthEast')
legend boxoff
titlefield='Basal Ice Melt';
filetag='BasalIceMelt';
title([titletag,' ',stime,' ',timetype,' ',titlefield])
ridgepack_fprint('png',[projtag,'_',titletag,'_',...
                        stime,'_',timetype,'_',filetag],1,2)


clf


ridgepack_polarm('antarctic','noland','grid','label');
geoshow(SCP,'FaceColor',0.95*[1 1 1],'EdgeColor',0.5*[1 1 1])
mask=find(nca.(iceareacell).data>0.001);
field=['time',timetype,'_avg_lateralIceMelt'];
nc=ridgepack_reduce(ridgepack_clone(datafile,...
                    {field},timeslice),{'time'});
ridgepack_e3smcolors(nc,field,ncvert,mask,...
                   [0 0.001 0.01 0.025 0.05 0.075 0.1:0.1:1.1]*10^-9);
h=geoshow(STC,'Color','m','LineWidth',1)
legend(h,'Extent','Location','SouthEast')
legend boxoff
titlefield='Lateral Ice Melt';
filetag='Lateral_Ice_Melt';
title([titletag,' ',stime,' ',timetype,' ',titlefield])
ridgepack_fprint('png',[projtag,'_',titletag,'_',...
                        stime,'_',timetype,'_',filetag],1,2)


clf

ridgepack_polarm('antarctic','noland','grid','label');
geoshow(SCP,'FaceColor',0.95*[1 1 1],'EdgeColor',0.5*[1 1 1])
mask=find(nca.(iceareacell).data>0.001);
field=['time',timetype,'_avg_snowMelt'];
nc=ridgepack_reduce(ridgepack_clone(datafile,...
                    {field},timeslice),{'time'});
ridgepack_e3smcolors(nc,field,ncvert,mask,...
                    [0 0.5 1:1:13]*10^-9);
h=geoshow(STC,'Color','m','LineWidth',1)
legend(h,'Extent','Location','SouthEast')
legend boxoff
titlefield='Snow Melt';
filetag='Snow_Melt';
title([titletag,' ',stime,' ',timetype,' ',titlefield])
ridgepack_fprint('png',[projtag,'_',titletag,'_',...
                        stime,'_',timetype,'_',filetag],1,2)


clf

ridgepack_polarm('antarctic','noland','grid','label');
geoshow(SCP,'FaceColor',0.95*[1 1 1],'EdgeColor',0.5*[1 1 1])
mask=find(nca.(iceareacell).data>0.001);
field=['time',timetype,'_avg_surfaceIceMelt'];
nc=ridgepack_reduce(ridgepack_clone(datafile,...
                    {field},timeslice),{'time'});
ridgepack_e3smcolors(nc,field,ncvert,mask,...
                    [0 0.001 0.01 0.1:0.1:1.1]*10^-9);
h=geoshow(STC,'Color','m','LineWidth',1)
legend(h,'Extent','Location','SouthEast')
legend boxoff
titlefield='Surface Ice Melt';
filetag='Surface_Ice_Melt';
title([titletag,' ',stime,' ',timetype,' ',titlefield])
ridgepack_fprint('png',[projtag,'_',titletag,'_',...
                        stime,'_',timetype,'_',filetag],1,2)


clf

ridgepack_polarm('antarctic','noland','grid','label');
geoshow(SCP,'FaceColor',0.95*[1 1 1],'EdgeColor',0.5*[1 1 1])
mask=find(nca.(iceareacell).data>0.001);
field=['time',timetype,'_avg_frazilFormation'];
nc=ridgepack_reduce(ridgepack_clone(datafile,...
                    {field},timeslice),{'time'});
ridgepack_e3smcolors(nc,field,ncvert,mask,...
                    [0:0.1:1.5]*10^-8);
h=geoshow(STC,'Color','m','LineWidth',1)
legend(h,'Extent','Location','SouthEast')
legend boxoff
titlefield='Frazil Formation';
filetag='Frazil_Formation';
title([titletag,' ',stime,' ',timetype,' ',titlefield])
ridgepack_fprint('png',[projtag,'_',titletag,'_',...
                        stime,'_',timetype,'_',filetag],1,2)


clf

ridgepack_polarm('antarctic','noland','grid','label');
geoshow(SCP,'FaceColor',0.95*[1 1 1],'EdgeColor',0.5*[1 1 1])
mask=find(nca.(iceareacell).data>=0.0);
field=['time',timetype,'_avg_seaSurfaceSalinity'];
nc=ridgepack_reduce(ridgepack_clone(datafile,...
                    {field},timeslice),{'time'});
ridgepack_e3smcolors(nc,field,ncvert,mask,...
                    [32.3:0.1:34.4]);
h=geoshow(STC,'Color','m','LineWidth',1)
legend(h,'Extent','Location','SouthEast')
legend boxoff
titlefield='Sea Surface Salinity';
filetag='Sea_Surface_Salinity';
title([titletag,' ',stime,' ',timetype,' ',titlefield])
ridgepack_fprint('png',[projtag,'_',titletag,'_',...
                        stime,'_',timetype,'_',filetag],1,2)


clf

ridgepack_polarm('antarctic','noland','grid','label');
geoshow(SCP,'FaceColor',0.95*[1 1 1],'EdgeColor',0.5*[1 1 1])
mask=find(nca.(iceareacell).data>=0.0);
field=['time',timetype,'_avg_seaSurfaceTemperature'];
nc=ridgepack_reduce(ridgepack_clone(datafile,...
                    {field},timeslice),{'time'});
ridgepack_e3smcolors(nc,field,ncvert,mask,...
                     [-1.8 -1.5:0.5:0 1:1:11]);
h=geoshow(STC,'Color','m','LineWidth',1)
legend(h,'Extent','Location','SouthEast')
legend boxoff
titlefield='Sea Surface Temperature';
filetag='Sea_Surface_Temperature';
title([titletag,' ',stime,' ',timetype,' ',titlefield])
ridgepack_fprint('png',[projtag,'_',titletag,'_',...
                        stime,'_',timetype,'_',filetag],1,2)


clf

ridgepack_polarm('antarctic','noland','grid','label');
geoshow(SCP,'FaceColor',0.95*[1 1 1],'EdgeColor',0.5*[1 1 1])
mask=find(nca.(iceareacell).data>0.001);
field=['time',timetype,'_avg_surfaceTemperatureCell'];
nc=ridgepack_reduce(ridgepack_clone(datafile,...
                    {field},timeslice),{'time'});
ridgepack_e3smcolors(nc,field,ncvert,mask,...
                    [-30:2:0]);
h=geoshow(STC,'Color','m','LineWidth',1)
legend(h,'Extent','Location','SouthEast')
legend boxoff
titlefield='Surface Temperature Cell';
filetag='Surface_Temperature_Cell';
title([titletag,' ',stime,' ',timetype,' ',titlefield])
ridgepack_fprint('png',[projtag,'_',titletag,'_',...
                        stime,'_',timetype,'_',filetag],1,2)


clf

ridgepack_polarm('antarctic','noland','grid','label');
geoshow(SCP,'FaceColor',0.95*[1 1 1],'EdgeColor',0.5*[1 1 1])
mask=find(nca.(iceareacell).data>0.001);
field=['time',timetype,'_avg_oceanHeatFlux'];
nc=ridgepack_reduce(ridgepack_clone(datafile,...
                    {field},timeslice),{'time'});
ridgepack_e3smcolors(nc,field,ncvert,mask,...
                     [-140:10:-10 -5 -2 -1 0],'linear',0);
h=geoshow(STC,'Color','m','LineWidth',1)
legend(h,'Extent','Location','SouthEast')
legend boxoff
titlefield='Ocean Heat Flux';
filetag='Ocean_Heat_Flux';
title([titletag,' ',stime,' ',timetype,' ',titlefield])
ridgepack_fprint('png',[projtag,'_',titletag,'_',...
                        stime,'_',timetype,'_',filetag],1,2)


clf

ridgepack_polarm('antarctic','noland','grid','label');
geoshow(SCP,'FaceColor',0.95*[1 1 1],'EdgeColor',0.5*[1 1 1])
mask=find(nca.(iceareacell).data>0.001);
field=['time',timetype,'_avg_congelation'];
nc=ridgepack_reduce(ridgepack_clone(datafile,...
                    {field},timeslice),{'time'});
ridgepack_e3smcolors(nc,field,ncvert,mask,...
                    [0:0.1:1.1]*10^-7);
h=geoshow(STC,'Color','m','LineWidth',1)
legend(h,'Extent','Location','SouthEast')
legend boxoff
titlefield='Congelation';
filetag='Congelation';
title([titletag,' ',stime,' ',timetype,' ',titlefield])
ridgepack_fprint('png',[projtag,'_',titletag,'_',...
                        stime,'_',timetype,'_',filetag],1,2)

clf

ridgepack_polarm('antarctic','noland','grid','label');
geoshow(SCP,'FaceColor',0.95*[1 1 1],'EdgeColor',0.5*[1 1 1])
mask=find(nca.(iceareacell).data>0.15);
field=['time',timetype,'_avg_iceAreaCell'];
nc=ridgepack_reduce(ridgepack_clone(datafile,...
                    {field},timeslice),{'time'});
ridgepack_e3smcolors(nc,field,ncvert,mask,...
                     [0:0.05:1]);
h=geoshow(STC,'Color','m','LineWidth',1)
legend(h,'Extent','Location','SouthEast')
legend boxoff
titlefield='Ice Area';
filetag='Ice_Area';
title([titletag,' ',stime,' ',timetype,' ',titlefield])
ridgepack_fprint('png',[projtag,'_',titletag,'_',...
                        stime,'_',timetype,'_',filetag],1,2)


clf

ridgepack_polarm('antarctic','noland','grid','label');
geoshow(SCP,'FaceColor',0.95*[1 1 1],'EdgeColor',0.5*[1 1 1])
mask=find(nca.(iceareacell).data>0.15);
field=['time',timetype,'_avg_iceVolumeCell'];
nc=ridgepack_reduce(ridgepack_clone(datafile,...
                    {field},timeslice),{'time'});
ridgepack_e3smcolors(nc,field,ncvert,mask,...
                     [0:0.25:3.25]);
h=geoshow(STC,'Color','m','LineWidth',1)
legend(h,'Extent','Location','SouthEast')
legend boxoff
titlefield='Ice Volume';
filetag='Ice_Volume';
title([titletag,' ',stime,' ',timetype,' ',titlefield])
ridgepack_fprint('png',[projtag,'_',titletag,'_',...
                        stime,'_',timetype,'_',filetag],1,2)


