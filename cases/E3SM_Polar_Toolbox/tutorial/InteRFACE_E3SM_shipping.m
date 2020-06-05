% InteRFACE_E3SM_shipping: Tutorial Script
% 
% This script demonstrates basic functionality for shipping.
%
% InteRFACE Tutorial, April 2020
% Andrew Roberts, LANL, afroberts@lanl.gov

clf
clear

% CHANGE THESE THINGS TO MAKE THIS WORK FOR YOU: 

% plot location
plotloc='/Users/afroberts/work/tutorial';

% grid data location
gridloc='/Users/afroberts/data/MODEL/E3SM/DECK/grid';
gridfile='E3SM_LR_V1_grid.nc';

% data file location
dataloc='/Users/afroberts/data/MODEL/E3SM/DECK/monthly/h1/archive/ice/hist';
datafile='mpascice.hist.am.timeSeriesStatsMonthly.2000-09-01.nc';

% EVERYTHING BELOW THIS LINE SHOULD JUST WORK

% get the mesh data
cd(gridloc)
ncvert=ridgepack_clone(gridfile,{'latVertex','lonVertex',...
                                'verticesOnCell','indexToCellID',...
                                'nEdgesOnCell','edgesOnCell',...
                                'cellsOnEdge','cellsOnVertex',...
                                'edgesOnVertex'});

nccell=ridgepack_clone(gridfile,{'latCell','lonCell'});

ncedge=ridgepack_clone(gridfile,{'latEdge','lonEdge'});

ncvert.latCell=nccell.latitude;
ncvert.lonCell=nccell.longitude;
ncvert.latEdge=ncedge.latitude;
ncvert.lonEdge=ncedge.longitude;

% grab a large chunk of different sea ice data, for demonstration
cd(dataloc)
fields={'timeMonthly_avg_iceAreaCell',...
        'timeMonthly_avg_iceVolumeCell'};

% ingest the data, and remove the time dimension immediately
nc=ridgepack_reduce(ridgepack_clone(datafile,fields),{'time'});

% Now, extract just the grid cells with this criteria
[nccoast,SCP,SCL,STP,STC]=ridgepack_e3smseasaw(ncvert,nc,...
                             'timeMonthly_avg_iceAreaCell',0.15);

% Reduce down the mesh

% Generate Great Circle Across the Arctic

% Bering Strait
lat1=67;
lon1=-169;

% North Sea
lat2=62;
lon2=3;

[dist,angl,phi,tracklat,tracklon,tracklen]=...
          ridgepack_greatcircle(lat1,lon1,lat2,lon2);

%globeplot=false;
globeplot=true;

if globeplot

 centerlat=70;
 centerlon=-45;
 horizon=90;

 ridgepack_satview(centerlat,centerlon,horizon)

 % Generate coast
 [x,y,z,phi,theta]=ridgepack_satfwd(nccoast.latitude.data,...
                                    nccoast.longitude.data,...
                                    centerlat,centerlon,horizon);

 plot3(x,y,z,'Color',0.5*[1 1 1],'LineWidth',0.75)

 [x,y,z,phi,theta]=ridgepack_satfwd(tracklat,tracklon,...
                            centerlat,centerlon,horizon);

 plot3(x,y,z,'Color','m','LineWidth',1.0)

 cd(plotloc)
 ridgepack_fprint('png',[mfilename,'_globe'],1,2)

else

 % create an Arctic base map
 ridgepack_polarm('shipping','noland','grid','label');

 % Plot sea ice thickness, masked at 15% concentration
 mask=find(nc.timeMonthly_avg_iceAreaCell.data>0.15);
 ridgepack_e3smcolors(nc,'timeMonthly_avg_iceVolumeCell',...
                     ncvert,mask,[0:0.25:3.25]);

 % Overlay extent line
 h=geoshow(STC,'Color',[0.88 0.3 0],'LineWidth',1.0);
 
 % Generate March Shipping Criteria of areas with less than 
 % 50% concentration.  All cells with this criteria will have 
 % value of 1, all others a value of zero.
 mask=find(nc.timeMonthly_avg_iceAreaCell.data<0.5);
 SML=ridgepack_e3smeshs(ncvert,mask,[0.88 0.3 0])

 % Plot the great circle ship path
 [x,y]= mfwdtran(gcm,tracklat,tracklon);
 hc=plot(x,y,'Color','m','LineWidth',1.0)

 % Add in the coast
 geoshow(SCL,'Color',0.5*[1 1 1])

 % Add Legend
 legend(hc,'Great Circle','Location','NorthWest')
 legend box off

 % title too
 title('Navigable Regions September 2000 E3SM DECK h1')

 % print it
 cd(plotloc)
 ridgepack_fprint('png',mfilename,1,2)

 % create a shape for the great circle
 SHIP=geoshape(tracklat,tracklon,...
             'Shipping','GreatCircle',...
             'Geometry','line');

 % save acessible mesh
 kmlwrite('Mesh_Shipping',SML,...
         'Color',[0.88 0.3 0],'LineWidth',0.5, ...
         'Description','Navigable Region','Name','MPAS Sea Ice');

 kmlwrite('Great_Circle',SHIP,'Color','w','LineWidth',0.5, ...
          'Description','Great Circle','Name','MPAS Sea Ice');

end

