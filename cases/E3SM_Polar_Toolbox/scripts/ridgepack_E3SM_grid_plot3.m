
clf
clear

tesselation=true;
%tesselation=false;

plotchoice=1;

for setting=plotchoice

 clf

 if setting==1

  centlat=-40; % degrees north
  centlon=-180; % degrees east
  horizon=5; % degrees of satellite horizon (0-90)
  altitude=1; % Mean Earth radius multiple
  titled='Pacific';

 end

 % plot location
 plotloc='/Users/afroberts/work';

 % grid location
 gridloc=['/Users/afroberts/data/MODEL/E3SM/v3/v3.LR.control/grid'];
 gridfile='20231209.v3.LR.piControl-spinup.chrysalis.mpaso.rst.1000-01-01_00000.nc';
 coastname='V3Icos'; % grid name

 % obtain grid information
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

 ncvert.latCell=nccell.latitude;
 ncvert.lonCell=nccell.longitude;
 ncvert.latEdge=ncedge.latitude;
 ncvert.lonEdge=ncedge.longitude;
 ncvert.latVertex=ncvert.latitude;
 ncvert.lonVertex=ncvert.longitude;
 
 nccell.latVertex=ncvert.latitude;
 nccell.lonVertex=ncvert.longitude;
 nccell.latEdge=ncedge.latitude;
 nccell.lonEdge=ncedge.longitude;
 nccell.latCell=nccell.latitude;
 nccell.lonCell=nccell.longitude;

 cd(plotloc)

 ridgepack_satview(centlat,centlon,horizon)

 if tesselation

  ridgepack_e3smsatmeshspecial(ncvert,centlat,centlon,horizon,altitude);

 else

  ridgepack_e3smsatmeshv(nccell,centlat,centlon,horizon,altitude);

 end

 %ridgepack_sathorizon(centlat,centlon,horizon,...
 %                       satlat,satlon,sathorizon,[0.83 0.5 0])

 title([titled,' ',coastname],'FontWeight','normal')

 if tesselation

  ridgepack_fprint('png',[coastname,'_tesselation_3_',num2str(setting)],1,2)

 else

  ridgepack_fprint('png',[coastname,'_triangulation_3_',num2str(setting)],1,2)

 end

end % setting












