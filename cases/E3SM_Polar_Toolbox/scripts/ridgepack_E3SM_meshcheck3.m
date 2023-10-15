
close all
clear

zoomedareas=true;
%zoomedareas=false;

%grids=[7 10];
%grids=[7 9];
grids=[7];
%grids=[8];
 
%maintitle='MPAS E3SM V3 Mesh';
maintitle='';

fileg{1}.name='WC14L64';
fileg{2}.name='WC14r03';
fileg{3}.name='DECK';
fileg{4}.name='ECwISC30to60E1r02';
fileg{5}.name='EC30to60E2r2';
fileg{6}.name='oQU480';
fileg{7}.name='ECwISC30to60E3r2';
fileg{8}.name='IcoswISC30E3r2';
fileg{9}.name='SOwISC12to60E3r1';
fileg{10}.name='RRSwISC6to18E3r1';
  
% plot location
plotloc='/Users/afroberts/work';
 
gridl=length(grids);
if gridl>2
 error('can only compare two different meshes')
end
  
for gridchoice=grids
 
   gridx=find(gridchoice==grids);
 
   % grid location
   if strcmp(char(fileg{gridchoice}.name),'DECK')
    gridloc=['/Users/afroberts/data/MODEL/E3SM/EC_60_30_Old/grid'];
    gridfile='init.nc';
   elseif strcmp(char(fileg{gridchoice}.name),'WC14L64')
    gridloc=['/Users/afroberts/data/MODEL/E3SM/WC14-64L/grid'];
    gridfile='initial_state.nc';
   elseif strcmp(char(fileg{gridchoice}.name),'WC14r03')
    gridloc=['/Users/afroberts/data/MODEL/E3SM/WC14/r03'];
    gridfile='initial_state.nc';
   elseif strcmp(char(fileg{gridchoice}.name),'ECwISC30to60E1r02')
    gridloc=['/Users/afroberts/data/MODEL/E3SM/ECwISC30to60E1r02'];
    gridfile='ocean.ECwISC30to60E1r02.200408.nc';
   elseif strcmp(char(fileg{gridchoice}.name),'EC30to60E2r2')
    gridloc=['/Users/afroberts/data/MODEL/E3SM/PIControlSI/grid'];
    gridfile='mpaso.rst.0002-01-01_00000.nc';
   elseif strcmp(char(fileg{gridchoice}.name),'oQU480')
    gridloc=['/Users/afroberts/work'];
    gridfile='mpassi.rst.0001-04-01_00000.nc';
   elseif strcmp(char(fileg{gridchoice}.name),'ECwISC30to60E3r2')
    gridloc=['/Users/afroberts/data/MODEL/E3SM/v3/ecwisc30to60e3r2'];
    gridfile='mpaso.ECwISC30to60E3r2.20230831.nc';
   elseif strcmp(char(fileg{gridchoice}.name),'IcoswISC30E3r2')
    gridloc=['/Users/afroberts/data/MODEL/E3SM/v3/IcoswISC30E3r2'];
    gridfile='mpaso.IcoswISC30E3r2.20230901.nc';
   elseif strcmp(char(fileg{gridchoice}.name),'SOwISC12to60E3r1')
    gridloc=['/Users/afroberts/data/MODEL/E3SM/v3/SOwISC12to60E3r1'];
    gridfile='mpaso.SOwISC12to60E3r1.20230901.nc';
   elseif strcmp(char(fileg{gridchoice}.name),'RRSwISC6to18E3r1')
    gridloc=['/Users/afroberts/data/MODEL/E3SM/v3/RRSwISC6to18E3r1'];
    gridfile='mpaso.RRSwISC6to18E3r1.20230902.nc';
   end
  
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
  
   % read in coast or else generate coast and write it out
   coastname=[char(fileg{gridchoice}.name),'_Coast.nc'];
   x=dir(coastname);
   if isempty(x)
    nccoast=ridgepack_e3smcoastm(ncvert);
    ridgepack_write(nccoast,coastname)
   else
    nccoast=ridgepack_clone(coastname);
   end
  
   % create isobaths
   bathname=[char(fileg{gridchoice}.name),'_30mIsobath.nc'];
   greename=[char(fileg{gridchoice}.name),'_50mIsobath.nc'];
   deepname=[char(fileg{gridchoice}.name),'_500mIsobath.nc'];
   x=dir(greename);
   if isempty(x) 
    ncisobath30=ridgepack_e3smseasaw(ncvert,ncvert,'bottomDepth',30);
    ridgepack_write(ncisobath30,bathname)
    ncisobath50=ridgepack_e3smseasaw(ncvert,ncvert,'bottomDepth',50);
    ridgepack_write(ncisobath50,bathname)
    ncisobath500=ridgepack_e3smseasaw(ncvert,ncvert,'bottomDepth',500);
    ridgepack_write(ncisobath500,bathname)
   else
    ncisobath30=ridgepack_clone(bathname);
    ncisobath50=ridgepack_clone(greename);
    ncisobath500=ridgepack_clone(deepname);
   end
  
end % gridchoice
   
% move to plot location
cd(plotloc)
  
title([char(fileg{gridchoice}.name)],'FontWeight','normal')


