function ridgepack_E3SM_resolution_check(centlat,centlon,...
                                         horizon,coastname,cont)

close all

% set constants
altitude=1; % Mean Earth radius multiple
cgrid=true; % plot c-grid coastline

% define contours
ref=0;

% grid location
if strcmp(coastname,'ARRM')
 %gridloc='/Users/afroberts/SIhMSatArray/E3SM/ARRM/grid';
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
 titlename='Greenland';
elseif strcmp(coastname,'CONUS_modified')
 gridloc='/Users/afroberts/SIhMSatArray/E3SM/Modified_CONUS/grid';
 gridfile='initial_state.nc';
 titlename='Modified CONUS';
elseif strcmp(coastname,'EC_60_30')
 gridloc='/Users/afroberts/SIhMSatArray/E3SM/EC_60_30/grid';
 gridfile='initial_state.nc';
 titlename='EC60to30';
elseif strcmp(coastname,'EC_60_30_Degraded')
 gridloc='/Users/afroberts/SIhMSatArray/E3SM/EC_60_30_Degraded/grid';
 gridfile='initial_state.nc';
 titlename='Degraded EC60to30';
elseif strcmp(coastname,'EC_60_30_Old')
 gridloc='/Users/afroberts/SIhMSatArray/E3SM/EC_60_30_Old/grid';
 gridfile='init.nc';
 titlename='Old EC60to30';
elseif strcmp(coastname,'DECK')
 gridloc='/Users/afroberts/data/E3SM/DECK/grid';
 gridfile='oEC60to30v3_60layer.restartFrom_anvil0926.171101.nc';
 titlename='DECK';
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

% plot cell resolution
ridgepack_multiplot(2,2,1,1,'a')
ridgepack_satview(centlat,centlon,90,1,0)
ridgepack_psatcole3sm(nccell,'areaCell',ncvert,cont,ref,...
                      centlat,centlon,90,1.001*altitude);
ridgepack_colorbar(cont,'km','linear','vertical',ref);
ridgepack_sathorizon(centlat,centlon,90,...
                     centlat,centlon,horizon,[0 0.8 0]);
ridgepack_psatcoaste3sm(ncvert,cgrid,coastname,...
                        centlat,centlon,90);
title([titlename,' $\sqrt{\textrm(Cell Area)}$'],'FontSize',10)

% plot mesh
ridgepack_multiplot(2,2,1,2,'b')
ridgepack_satview(centlat,centlon,horizon,1,0);
ridgepack_psatmeshe3sm(nccell,'areaCell',ncvert,...
           cont,ref,centlat,centlon,horizon,altitude);
ridgepack_psatcoaste3sm(ncvert,cgrid,coastname,...
                        centlat,centlon,horizon);
title('Analysis Region','FontSize',10)

% get resolution statistics for the grid cells
ridgepack_multiplot(2,2,2,1,'c')
[cells]=ridgepack_psatmeshe3sm(nccell,'areaCell',ncvert,...
                   cont,ref,centlat,centlon,horizon,altitude);

dc=[];
for i=1:length(cells)
 idx=ncvert.maxEdges.data(1:ncvert.nEdgesOnCell.data(cells(i)));
 edgeidx=ncvert.edgesOnCell.data(idx,cells(i));
 dc=[dc; ncvert.dcEdge.data(edgeidx)/1000]; % convert to km
end

[X,Y]=hist(dc,100); % obtain histogram
X=X(:)./sum(X(:));
stairs(Y,X)
xlim([min(Y) max(Y)])
xlabel('Cell-to-Cell Resolution (km)','FontSize',10)
if centlat<0
 latunits='S';
else
 latunits='N';
end
if centlon<0
 lonunits='W';
else
 lonunits='E';
end
title(['Edges on cells within ',num2str(horizon),...
        '$^{\circ}$'],'FontSize',10)

% get resolution statistics for the grid cells
ridgepack_multiplot(2,2,2,2,'d')
[cells]=ridgepack_psatmeshe3sm(nccell,'areaCell',ncvert,...
                   cont,ref,centlat,centlon,horizon,altitude);

diffc=NaN*zeros(size(cells));

for i=1:length(cells)

 idx=ncvert.maxEdges.data(1:ncvert.nEdgesOnCell.data(cells(i)));

 edgeidx=ncvert.edgesOnCell.data(idx,cells(i));

 diffc(i)=100*(max(ncvert.dcEdge.data(edgeidx))-...
               min(ncvert.dcEdge.data(edgeidx)))./...
               min(ncvert.dcEdge.data(edgeidx));

end

[X,Y]=hist(diffc,50); % obtain histogram
X=X(:)./sum(X(:));
Y=Y(:); % convert to km
stairs(Y,X)
xlim([0 max(Y)])
xlabel('Cell-to-Cell Resolution Change (\%)','FontSize',10)
if centlat<0
 latunits='S';
else
 latunits='N';
end
if centlon<0
 lonunits='W';
else
 lonunits='E';
end
title(['Cells within ',num2str(horizon),...
        '$^{\circ}$ of ',...
        num2str(abs(centlat)),latunits,' ',...
        num2str(abs(centlon)),lonunits],'FontSize',10)

ridgepack_multialign(gcf)

cd ~/Work
plotname=[coastname,'_',num2str(horizon),'_',...
         num2str(abs(centlat)),latunits,'_',...
         num2str(abs(centlon)),lonunits,'_gridcheck'];

ridgepack_fprint('png',plotname,1,1)
ridgepack_fprint('epsc',plotname,1,1)


