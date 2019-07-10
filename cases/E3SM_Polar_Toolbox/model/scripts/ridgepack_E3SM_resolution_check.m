function ridgepack_E3SM_resolution_check(centlat,centlon,...
                                         horizon,coastname,cont)

close all

% Arctic Center
%centlat=50; % degrees north
%centlon=-20; % degrees east
%horizon=10; % degrees of satellite horizon (0-90)
%coastname='Greenland'; % grid name
%cont=[10:5:60];

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
elseif strcmp(coastname,'Greenland')
 gridloc='/Users/afroberts/data/E3SM/Greenland/grid';
 gridfile='greenland_grid.nc';
elseif strcmp(coastname,'CONUS')
 gridloc='/Users/afroberts/data/E3SM/CONUS/grid';
 gridfile='initial_state.nc';
elseif strcmp(coastname,'DECK')
 gridloc='/Users/afroberts/data/E3SM/DECK/grid';
 gridfile='oEC60to30v3_60layer.restartFrom_anvil0926.171101.nc';
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
title([coastname,' $\sqrt{\textrm(areaCell)}$'],'FontSize',10)

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
idx=ncvert.edgesOnCell.data(:,cells);
id=unique(sort(idx(idx(:)>0)));
dc=ncvert.dcEdge.data(id);
[X,Y]=hist(dc,100); % obtain histogram
X=X(:)./sum(X(:));
Y=Y(:)./1000; % convert to km
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
 idx=find(ncvert.edgesOnCell.data(:,i)>0);  
 edgeidx=ncvert.edgesOnCell.data(idx,i);
 diffc(i)=1-(min(ncvert.dcEdge.data(edgeidx))./...
             max(ncvert.dcEdge.data(edgeidx)));
end

[X,Y]=hist(diffc,50); % obtain histogram
X=X(:)./sum(X(:));
Y=Y(:).*100; % convert to km
stairs(Y,X)
xlim([0 max(Y)])
xlabel('Cell-to-Cell Resolution Gradient (\%)','FontSize',10)
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


