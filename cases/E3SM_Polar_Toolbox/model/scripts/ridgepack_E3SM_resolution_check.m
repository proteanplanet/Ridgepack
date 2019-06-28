
clf
clear

% Arctic Center
centlat=79; % degrees north
centlon=-100; % degrees east
horizon=3; % degrees of satellite horizon (0-90)
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
                                 'cellsOnEdge','areaCell',...
                                 'weightsOnEdge'});

nccell=ridgepack_clone(gridfile,{'latCell','lonCell','areaCell'});

nccell.areaCell.data=sqrt(nccell.areaCell.data)/1000;

% plot cell resolution
ridgepack_multiplot(1,3,1,1,'a')
ridgepack_satview(centlat,centlon,90,1,1)
ridgepack_psatcole3sm(nccell,'areaCell',ncvert,cont,ref,...
                      centlat,centlon,90,altitude);

% plot mesh
ridgepack_multiplot(1,3,1,2,'b')
ridgepack_satview(centlat,centlon,horizon,1,1)
ridgepack_psatmeshe3sm(nccell,'areaCell',ncvert,...
                   cont,ref,centlat,centlon,horizon,altitude);

% plot coastal outline
ridgepack_psatcoaste3sm(ncvert,cgrid,coastname,...
                             centlat,centlon,horizon);

% get resolution statistics for the grid cells
ridgepack_multiplot(1,3,1,3,'c')
[cells]=ridgepack_psatmeshe3sm(nccell,'areaCell',ncvert,...
                   cont,ref,centlat,centlon,horizon,altitude);
idx=ncvert.edgesOnCell.data(:,cells);
id=unique(sort(idx(:)));
dc=ncvert.dcEdge.data(id);
[X,Y]=hist(dc,100); % obtain histogram
Y=Y/1000; % convert to km
stairs(Y,X)
xlabel('cell-to-cell resolution (km)')
if centlat<0
 latunits='$^{\circ}$S '
else
 latunits='$^{\circ}$N '
end
if centlon<0
 lonunits='$^{\circ}$W'
else
 lonunits='$^{\circ}$E'
end
title(['Edges within ',num2str(horizon),'$^{\circ}$ of ',...
        num2str(abs(centlat)),latunits,num2str(abs(centlon)),lonunits])

ridgepack_multialign(gcf)


cd ~/Work
ridgepack_fprint('png','test',1,1)


