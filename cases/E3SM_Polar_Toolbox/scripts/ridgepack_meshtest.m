
clf
clear

centlat=65; % degrees north
centlon=-40; % degrees east
horizon=8; % degrees of satellite horizon (0-90)
altitude=1; % Mean Earth radius multiple

% plot location
plotloc='/Users/afroberts/Documents';

% grid location
gridloc=['/Users/afroberts/data/MODEL/E3SM/WC14/r03'];
gridfile='initial_state.nc';

% obtain grid information
cd(gridloc)
ncvert=ridgepack_clone(gridfile,{'latVertex','lonVertex',...
                                'verticesOnCell','indexToCellID',...
                                'nEdgesOnCell','edgesOnCell',...
                                'cellsOnEdge','cellsOnVertex',...
                                'edgesOnVertex','bottomDepth'});

nccell=ridgepack_clone(gridfile,{'latCell','lonCell',...
                                'verticesOnCell','indexToCellID',...
                                'nEdgesOnCell','edgesOnCell',...
                                'cellsOnEdge','cellsOnVertex',...
                                'edgesOnVertex','bottomDepth'});

ncedge=ridgepack_clone(gridfile,{'latEdge','lonEdge'});

ncvert.latCell=nccell.latitude;
ncvert.lonCell=nccell.longitude;
ncvert.latEdge=ncedge.latitude;
ncvert.lonEdge=ncedge.longitude;
ncvert.latVertex=ncvert.latitude;
ncvert.lonVertex=ncvert.longitude;

% move to plot location
%mkdir(plotloc)
cd(plotloc)

ridgepack_satview(centlat,centlon,horizon)

ridgepack_e3smsatmeshs(ncvert,centlat,centlon,horizon,altitude);

ridgepack_fprint('png','E3SM_v2_mesh_plot_test',1,2)

