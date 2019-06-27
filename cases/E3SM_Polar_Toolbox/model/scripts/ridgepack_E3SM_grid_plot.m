
clf
clear

% Arctic Center
centlat=79; % degrees north
centlon=-100; % degrees east
horizon=90; % degrees of satellite horizon (0-90)
altitude=1; % Mean Earth radius multiple

% define contours
cont=[5 10:1:31];
ref=20;

% grid location
%gridloc='/Users/afroberts/SIhMSatArray/E3SM/ARRM/grid';
gridloc='/Users/afroberts/data/E3SM/ARRM/grid';

% obtain grid information
cd(gridloc)
gridfile='ocean.ARRM60to10.180715.nc';
ncvert=ridgepack_clone(gridfile,{'latVertex','lonVertex',...
                                 'verticesOnCell','indexToCellID',...
                                 'nEdgesOnCell','areaCell'});

nccell=ridgepack_clone(gridfile,{'latCell','lonCell',...
                                 'verticesOnCell','indexToCellID',...
                                 'nEdgesOnCell','areaCell'});

nccell.areaCell.data=sqrt(nccell.areaCell.data)/1000;

ridgepack_satview(centlat,centlon,horizon,1,2)

%ridgepack_psatcole3sm(nccell,'areaCell',ncvert,nccell,cont,ref,...
%                      centlat,centlon,horizon,altitude);

ridgepack_psatmeshe3sm(nccell,'areaCell',ncvert,nccell,cont,ref,...
                      centlat,centlon,horizon,altitude);

% obtain coastal outline
cd(gridloc)
landm=shaperead('E3SM_ARRM_V1_C_grid_Coast.shp','UseGeoCoords',true);


