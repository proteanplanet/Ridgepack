
clear
clf

% plot location
centlat=40;
centlon=18;
horizon=10.0;
altitude=1;
name='Italy';

% add info including tag 'Critical_Passage' or 'Critical_Land_Blockage'
track{1}.name='Sicily';
track{1}.tag='Critical_Land_Blockage';

track{2}.name='Calabria';
track{2}.tag='Critical_Land_Blockage';

track{3}.name='Salento';
track{3}.tag='Critical_Land_Blockage';

% added grid
gridchoice=2;

fileg{1}.name='WC12r01';
fileg{1}.outname='WC12';
fileg{1}.title=' WC 12-60~km mesh';
fileg{2}.name='WC14r03';
fileg{2}.outname='WC14r03';
fileg{2}.title=' WC 14-60~km mesh r03';
fileg{3}.name='DECK';
fileg{3}.outname='DECK';
fileg{3}.title=' DECK 30-60~km standard mesh';

% grid location
if strcmp(char(fileg{gridchoice}.name),'DECK')
 gridloc=['/Users/afroberts/data/MODEL/E3SM/',...
          char(fileg{gridchoice}.name),'/grid'];
 gridfile='E3SM_LR_V1_grid.nc';
elseif strcmp(char(fileg{gridchoice}.name),'WC12')
 gridloc=['/Users/afroberts/data/MODEL/E3SM/WC12/',...
          char(fileg{gridchoice}.name)];
 gridfile='initial_state.nc';
elseif strcmp(char(fileg{gridchoice}.name),'WC14r03')
 gridloc=['/Users/afroberts/data/MODEL/E3SM/WC14/r03'];
 gridfile='initial_state.nc';
end

gridlochr=['/Users/afroberts/data/MODEL/E3SM/highres/grid'];
gridfilehr='E3SM_hr_grid.nc';

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
coastname=[char(fileg{gridchoice}.outname),'_Coast.nc'];
x=dir(coastname);
if isempty(x)
 nccoast=ridgepack_e3smcoastm(ncvert);
 ridgepack_write(nccoast,coastname)
else
 nccoast=ridgepack_clone(coastname);
end

% create 20m isobath
bathname=[char(fileg{gridchoice}.outname),'_20mIsobath.nc'];
x=dir(bathname);
if isempty(x)
  ncisobath20=ridgepack_e3smseasaw(ncvert,ncvert,'bottomDepth',20);
  ridgepack_write(ncisobath20,bathname)
else
  ncisobath20=ridgepack_clone(bathname);
end

% invert bathymetry
ncvert.bottomDepth.data=-ncvert.bottomDepth.data;

% load high-resolution coast
cd(gridlochr)
coastnamehr='E3SM_HR_V1_Coast';
nccoasthr=ridgepack_clone(coastnamehr);

shiploc='/Users/afroberts/data/Transects';
cd(shiploc)

ocean=false;
land=false;

for k=1:length(track)

 if strcmp(char(track{k}.tag),'Critical_Land_Blockage')
  land=true;
  thisland=true;
 elseif strcmp(char(track{k}.tag),'Critical_Passage')
  ocean=true;
 else
  thisland=false;
 end

 kmlStruct=kml2struct([char(track{k}.name),'.kml']);
 coords{k}.lats=kmlStruct.Lat;
 coords{k}.lons=kmlStruct.Lon;

 % write out GeoJSON
 filename=[char(track{k}.name),'.geojson'];
 fileID = fopen(filename,'w');
 if thisland
  fprintf(fileID,['{"type":"FeatureCollection","features":[{"type":"Feature","properties":{"name":"',char(track{k}.name),'","tags":"',char(track{k}.tag),'","component":"ocean","author":"Andrew Roberts","object":"transect","height":"100.0"},"geometry":{"type":"LineString","coordinates":[']);
 else
  fprintf(fileID,['{"type":"FeatureCollection","features":[{"type":"Feature","properties":{"name":"',char(track{k}.name),'","tags":"',char(track{k}.tag),'","component":"ocean","author":"Andrew Roberts","object":"transect"},"geometry":{"type":"LineString","coordinates":[']);
 end
 for j=1:length(coords{k}.lats)-1
  fprintf(fileID,['[',num2str(coords{k}.lons(j),'%2.10f'),',',...
                      num2str(coords{k}.lats(j),'%2.10f'),'],']);
 end
 fprintf(fileID,['[',num2str(coords{k}.lons(j),'%2.10f'),',',...
                       num2str(coords{k}.lats(j),'%2.10f'),']']);
 fprintf(fileID,']}}]}');
 fclose(fileID);

 [coords{k}.x,coords{k}.y,coords{k}.z]=...
                          ridgepack_satfwd(coords{k}.lats,...
                                           coords{k}.lons,...
                          centlat,centlon,horizon,1.001*altitude);

end

ridgepack_multiplot(1,2,1,1)
ridgepack_satview(centlat,centlon,horizon)
ridgepack_e3smsatmeshs(ncvert,centlat,centlon,horizon,altitude);

if land & ocean
 for k=1:length(track)
  hland=plot3(coords{k}.x,coords{k}.y,coords{k}.z,'r-');
  hocean=plot3(coords{k}.x,coords{k}.y,coords{k}.z,'m-');
 end
elseif land
 for k=1:length(track)
  hland=plot3(coords{k}.x,coords{k}.y,coords{k}.z,'r-');
 end
elseif ocean
 for k=1:length(track)
  hocean=plot3(coords{k}.x,coords{k}.y,coords{k}.z,'m-');
 end
end

% add 20m isobath
hk=ridgepack_e3smsatthreshold(ncisobath20,centlat,centlon,horizon,...
                                 [0 0 1]);
set(hk,'LineWidth',0.5)

% add in coastline
ridgepack_e3smsatcoast(nccoast,centlat,centlon,horizon)

% add in high resolution coastline
h=ridgepack_e3smsatcoast(nccoasthr,centlat,centlon,horizon,[0 0.8 0])

ridgepack_multiplot(1,2,1,2)

ridgepack_satview(centlat,centlon,horizon,1,1)

% reverse colorbar
cont=[-5000:1000:-2000 -1500 -1000:250:-250 -100 -50:10:10];
colbarcont{1}='\downarrow';
for i=2:length(cont)-1;
  colbarcont{i}=num2str(-cont(i));
end
colbarcont{length(cont)}='\uparrow';
cmap=colormap;

% render colors
ridgepack_e3smsatcol(ncvert,'bottomDepth',ncvert,cont,0,...
                        centlat,centlon,horizon,altitude,...
                        true,false,'linear','bluered');

% make top color grey
cmap=colormap;
cmap(end,:)=0.95*[1 1 1];
colormap(cmap)

% add colormap
ridgepack_colorbar(cont,'m','linear','vertical',0,colbarcont)
clear cmap colbarcont

% add 20m isobath
hb=ridgepack_e3smsatthreshold(ncisobath20,centlat,centlon,horizon,...
                                 [0.9290 0.6940 0.1250]);

% add coast
ridgepack_e3smsatcoast(nccoast,centlat,centlon,horizon)

% added transects
if land & ocean
 for k=1:length(track)
  hland=plot3(coords{k}.x,coords{k}.y,coords{k}.z,'r-');
  hocean=plot3(coords{k}.x,coords{k}.y,coords{k}.z,'m-');
 end
 ridgepack_multilegend([h hb hland hocean],...
   {'E3SM-HR V1 coastline','20 m isobath','New Enforced Land','New Enforced Passage'},'South')
elseif land 
 for k=1:length(track)
  hland=plot3(coords{k}.x,coords{k}.y,coords{k}.z,'r-');
 end
 ridgepack_multilegend([h hb hland],...
   {'E3SM-HR V1 coastline','20 m isobath','New Enforced Land'},'South')
elseif ocean 
 for k=1:length(track)
  hocean=plot3(coords{k}.x,coords{k}.y,coords{k}.z,'m-');
 end
 ridgepack_multilegend([h hb hocean],...
   {'E3SM-HR V1 coastline','20 m isobath','New Enforced Passage'},'South')
end



ridgepack_multialign(gcf,...
             [name,' ',char(fileg{gridchoice}.title)]);

ridgepack_fprint('png',[name,'_',char(fileg{gridchoice}.name),'_fix'],1,2)


