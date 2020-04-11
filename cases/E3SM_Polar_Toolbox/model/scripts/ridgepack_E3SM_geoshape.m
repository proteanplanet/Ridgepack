
clf
clear

regenerate=true;
%regenerate=false;

debug=false;
%debug=true;

centlat=80; % degrees north
centlon=-150; % degrees east
horizon=90; % degrees of satellite horizon (0-90)

% location of grid file
gridloc='/Users/afroberts/data/MODEL/E3SM/DECK/grid';
gridfile='E3SM_LR_V1_grid.nc';
cgrid=true;

% location of sea ice data
dataloc='/Users/afroberts/data/MODEL/E3SM/DECK/monthly/h1/archive/ice/reduced';
datafile='mpascice.hist.am.timeSeriesStatsMonthly.1980-03-01.nc';
varc='timeMonthly_avg_iceAreaCell';
datatitle='Sea Ice Extent';
threshold=0.15;

% plot location
plotloc='/Users/afroberts/work';

% %%%%%%%%%%%%%%% CHANGE ABOVE THIS LINE %%%%%%%%%%%%%%%%%% %


% obtain grid information (vertices, cell centers)
if regenerate

 cd(gridloc)
 ncvert=ridgepack_clone(gridfile,{'latVertex','lonVertex','dcEdge',...
                                 'verticesOnCell','indexToCellID',...
                                 'nEdgesOnCell','edgesOnCell',...
                                 'cellsOnEdge'});

 nccell=ridgepack_clone(gridfile,{'latCell','lonCell'});

 % obtain data
 cd(dataloc)
 ncc=ridgepack_clone(datafile,varc);

 % put up satview
 %ridgepack_satview(centlat,centlon,horizon,1,2);

 % plot ice edge 
 if 1==0
  ridgepack_pmeshbounde3sm(ncvert,ncc,varc,threshold,...
                         centlat,centlon,horizon);
 else
  nc=ridgepack_pmeshbounde3sm(ncvert,ncc,varc,threshold);

  cd(plotloc)
  save('testdata','nc','ncvert')
 end

else

 cd(plotloc)
 load testdata

end


S=geoshape(nc.latitude.data,nc.longitude.data,...
            'MPAS_Threshold','Sea Ice Threshold');

kmlwrite('MPAS_Threshold',S,'Color','m','Name','Extent',...
         'Description','Sea Ice Threshold');

clf

ridgepack_satview(centlat,centlon,horizon,1,2);

[x,y,z,phi,theta]=ridgepack_satfwd(nc.latitude.data,...
                                   nc.longitude.data,...
                      centlat,centlon,2*horizon,1.0001);
plot3(x,y,z,'Color','m',...
       'LineWidth',0.5-sin(deg2rad(horizon))*0.4)

hold on

[x,y,z,phi,theta]=ridgepack_satfwd(nc.clatitude.data,...
                                   nc.clongitude.data,...
                      centlat,centlon,2*horizon,1.0001);
plot3(x,y,z,'Color','k',...
       'LineWidth',0.5-sin(deg2rad(horizon))*0.4)

hold on
hold on

 
% add title
%title([datatitle])

cd(plotloc)
ridgepack_fprint('png','Outfile_Sea_Ice_GeoTif',1,1)




