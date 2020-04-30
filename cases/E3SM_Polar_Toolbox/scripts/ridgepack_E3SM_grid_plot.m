
clf
clear

% Arctic Center
centlat=-69; % degrees north
centlon=-100; % degrees east
horizon=10; % degrees of satellite horizon (0-90)
altitude=1; % Mean Earth radius multiple
cgrid=true; % plot c-grid coastline
coastname='ARRM'; % grid name

% define contours
cont=[5 10:1:31];
ref=20;

% grid location
gridloc='/Users/afroberts/SIhMSatArray/E3SM/ARRM/grid';
%gridloc='/Users/afroberts/data/E3SM/ARRM/grid';

% obtain grid information
cd(gridloc)
gridfile='ocean.ARRM60to10.180715.nc';
ncvert=ridgepack_clone(gridfile,{'latVertex','lonVertex','dcEdge',...
                                 'verticesOnCell','indexToCellID',...
                                 'nEdgesOnCell','edgesOnCell',...
                                 'cellsOnEdge','areaCell'});

nccell=ridgepack_clone(gridfile,{'latCell','lonCell','areaCell'});

nccell.areaCell.data=sqrt(nccell.areaCell.data)/1000;

ridgepack_satview(centlat,centlon,horizon,1,2)

% plot cell resolution
%ridgepack_psatcole3sm(nccell,'areaCell',ncvert,nccell,cont,ref,...
%                      centlat,centlon,horizon,altitude);

% plot mesh
[cells]=ridgepack_psatmeshe3sm(nccell,'areaCell',ncvert,...
                   cont,ref,centlat,centlon,horizon,altitude);

% plot coastal outline
ridgepack_psatcoaste3sm(ncvert,cgrid,coastname,...
                             centlat,centlon,horizon);

cd ~/Work
ridgepack_fprint('png','test',1,1)


