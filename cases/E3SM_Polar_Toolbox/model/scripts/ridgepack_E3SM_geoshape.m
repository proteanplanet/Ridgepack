
clf
clear

regenerate=true;
%regenerate=false;

debug=false;
%debug=true;

centlat=90; % degrees north
centlon=0; % degrees east
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
 ncvert=ridgepack_clone(gridfile,{'latVertex','lonVertex',...
                                 'dcEdge',...
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
  ridgepack_e3smseasaw(ncvert,ncc,varc,threshold,...
                         centlat,centlon,horizon);
 else
  nc=ridgepack_e3smseasaw(ncvert,ncc,varc,threshold);

  cd(plotloc)
  save('testdata','nc','ncvert')
 end

else

 cd(plotloc)
 load testdata

end

S=geoshape(nc.latitude.data,nc.longitude.data,...
            'MPAS_Threshold','Sea Ice Threshold');

SP=geoshape(nc.latitude.data,nc.longitude.data,...
            'MPAS_Threshold','Sea Ice Threshold',...
            'Geometry','polygon');

ST=geoshape(nc.latitude.data,nc.longitude.data,...
            'MPAS','Sea Ice Threshold',...
            'Geometry','polygon');

SC=geoshape(nc.clatitude.data,nc.clongitude.data,...
            'MPAS','Sea Ice Threshold',...
            'Geometry','polygon');

kmlwrite('MPAS_Threshold',ST,...
         'Color',[1 1 1],...
         'Name','Extent',...
         'Description','Sea Ice Threshold');

kmlwrite('MPAS_Land',SC,...
         'FaceColor',0.9*[1 1 1],...
         'EdgeColor',0.5*[1 1 1],...
         'Name','Land',...
         'Description','Sea Ice Threshold');

clf

ridgepack_satview(centlat,centlon,horizon,1,2);

[x,y,z,phi,theta]=ridgepack_satfwd(nc.latitude.data,...
                                   nc.longitude.data,...
                      centlat,centlon,2*horizon,1.0001);
plot3(x,y,z,'Color','g',...
       'LineWidth',0.5)

hold on


[x,y,z,phi,theta]=ridgepack_satfwd(nc.xlatitude.data,...
                                   nc.xlongitude.data,...
                      centlat,centlon,2*horizon,1.0001);
plot3(x,y,z,'Color','m',...
       'LineWidth',0.5)

hold on

[x,y,z,phi,theta]=ridgepack_satfwd(nc.clatitude.data,...
                                   nc.clongitude.data,...
                      centlat,centlon,2*horizon,1.0001);
plot3(x,y,z,'Color','k',...
       'LineWidth',0.1)

hold on
 
% add title
%title([datatitle])

cd(plotloc)
ridgepack_fprint('png','Outfile_Sea_Ice_GeoTif',1,2)




