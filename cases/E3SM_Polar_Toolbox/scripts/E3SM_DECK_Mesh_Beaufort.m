
clf
clear 

withmesh=true;
%withmesh=false;

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

[nccoast,SCP]=ridgepack_e3smseasaw(ncvert)

ridgepack_polarm('centralarctic','noland','grid');
ridgepack_e3smeshv(nccell)
ridgepack_e3smeshs(ncvert)
plot([SCP.Latitude],[SCP.Longitude],'Color','k','LineWidth',0.75)
title('Arctic Ocean Delaunay Triangulation')
plotloc='/Users/afroberts/work';
cd(plotloc)
ridgepack_fprint('png',...
   ['Arctic_Delaunay_Triangulation'],1,2)

clf

ridgepack_polarm('lat',88,'noland');
ridgepack_e3smeshv(nccell)
ridgepack_e3smeshs(ncvert)
plot([SCP.Latitude],[SCP.Longitude],'Color','k','LineWidth',0.75)
title('North Pole Delaunay Triangulation')
plotloc='/Users/afroberts/work';
cd(plotloc)
ridgepack_fprint('png',...
   ['North_Pole_Delaunay_Triangulation'],1,2)


