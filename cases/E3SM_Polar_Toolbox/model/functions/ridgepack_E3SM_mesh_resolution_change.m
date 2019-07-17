function ridgepack_E3SM_mesh_resolution_change(coastname)

cont=[0:5:55];

close all

% set constants
altitude=1; % Mean Earth radius multiple
cgrid=true; % plot c-grid coastline

% define contours
ref=15;

% grid location
if strcmp(coastname,'ARRM')
 gridloc='/Users/afroberts/data/E3SM/ARRM/grid';
 gridfile='ocean.ARRM60to10.180715.nc';
 titlename='ARRM';
elseif strcmp(coastname,'Greenland')
 gridloc='/Users/afroberts/data/E3SM/Greenland/grid';
 gridfile='greenland_grid.nc';
 titlename='Greenland';
elseif strcmp(coastname,'CONUS')
 gridloc='/Users/afroberts/data/E3SM/CONUS/grid';
 gridfile='initial_state.nc';
 titlename='CONUS';
elseif strcmp(coastname,'CONUS_modified')
 gridloc='/Users/afroberts/data/E3SM/Modified_CONUS/grid';
 gridfile='initial_state.nc';
 titlename='Modified CONUS';
elseif strcmp(coastname,'EC_60_30')
 gridloc='/Users/afroberts/data/E3SM/EC_60_30/grid';
 gridfile='initial_state.nc';
 titlename='EC60to30';
elseif strcmp(coastname,'EC_60_30_Degraded')
 gridloc='/Users/afroberts/data/E3SM/EC_60_30_Degraded/grid';
 gridfile='initial_state.nc';
 titlename='Degraded EC60to30';
elseif strcmp(coastname,'EC_60_30_Old')
 gridloc='/Users/afroberts/data/E3SM/EC_60_30_Old/grid';
 gridfile='init.nc';
 titlename='Old EC60to30';
elseif strcmp(coastname,'DECK')
 gridloc='/Users/afroberts/data/E3SM/DECK/grid';
 gridfile='oEC60to30v3_60layer.restartFrom_anvil0926.171101.nc';
 titlename='DECK';
else
 error('coastname not recognized')
end

% obtain grid information
cd(gridloc)
ncvert=ridgepack_clone(gridfile,{'latVertex','lonVertex',...
                                 'dcEdge',...
                                 'verticesOnCell',...
                                 'indexToCellID',...
                                 'nEdgesOnCell',...
                                 'edgesOnCell',...
                                 'cellsOnEdge',...
                                 'areaCell',...
                                 'weightsOnEdge'});

nccell=ridgepack_clone(gridfile,{'latCell',...
                                 'lonCell',...
                                 'areaCell'});

nccell.areaCell.data=sqrt(nccell.areaCell.data)/1000;


cells=ncvert.nCells.data;

diffc=NaN*zeros(size(cells));

for i=1:length(cells)

 idx=ncvert.maxEdges.data(1:ncvert.nEdgesOnCell.data(cells(i)));

 edgeidx=ncvert.edgesOnCell.data(idx,cells(i));

 diffc(i)=100*(max(ncvert.dcEdge.data(edgeidx))-...
               min(ncvert.dcEdge.data(edgeidx)))./...
               min(ncvert.dcEdge.data(edgeidx));

end

nccell.diffc=nccell.areaCell;
nccell.diffc.units='percent';
nccell.diffc.long_name='Percent change in resolution across a cell';
nccell.diffc.data=diffc';

lat=[0 0 90 -90];
lon=[0 180 0 0];

kk=0;

for row=1:2
for col=1:2
 
 kk=kk+1;

 ridgepack_multiplot(2,2,row,col)

 ridgepack_satview(lat(kk),lon(kk),90,1,0)

 ridgepack_psatcole3sm(nccell,'diffc',ncvert,cont,ref,...
                      lat(kk),lon(kk),90,1.001*altitude);

 if kk==1
   ridgepack_colorbar(cont,'\%','linear','vertical',ref);
   ridgepack_multicb(gca)
 end

 ridgepack_psatcoaste3sm(ncvert,cgrid,coastname,...
                        lat(kk),lon(kk),90);

end
end

ridgepack_multialign(gcf,[titlename,...
            ' Intercell Maximum Resolution Change'])

cd ~/Work

plotname=[coastname,'_reschange'];

ridgepack_fprint('png',plotname,1,1)
ridgepack_fprint('epsc',plotname,1,1)

