
clf
clear

% set plot configuration
centlat=65; % degrees north
centlon=-40; % degrees east
horizon=8; % degrees of satellite horizon (0-90)
altitude=1; % Mean Earth radius multiple

% set plot location
plotloc='/Users/afroberts/Documents';

% set MPAS restart file location
gridloc=['/Users/afroberts/data/MODEL/E3SM/WC14/r03'];
gridfile='initial_state.nc';

% obtain grid information from the restart file
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

% add necessary information to cell array
ncvert.latCell=nccell.latitude;
ncvert.lonCell=nccell.longitude;
ncvert.latEdge=ncedge.latitude;
ncvert.lonEdge=ncedge.longitude;
ncvert.latVertex=ncvert.latitude;
ncvert.lonVertex=ncvert.longitude;

% move to plot location
%mkdir(plotloc)
cd(plotloc)

% plot satellite horizon
ridgepack_satview(centlat,centlon,horizon)

% plot MPAS mesh
ridgepack_e3smsatmeshs(ncvert,centlat,centlon,horizon,altitude);

% print the figure
ridgepack_fprint('png','E3SM_v2_mesh_plot_test',1,2)

