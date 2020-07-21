
clear

shiploc='/Users/afroberts/data/Transects';

clf

cd(shiploc)

kmlStruct=kml2struct('Foxe_Basin_Through_Flow.kml');

k=1
trackname{k}='Foxe_Basin_Through_Flow';
coords{k}.lats=kmlStruct.Lat;
coords{k}.lons=kmlStruct.Lon;

ridgepack_polarm('shipping')

for k=1:length(trackname)

 % write out GeoJSON
 filename=['E3SM_Foxe_Basin_Through_Flow.geojson'];
 fileID = fopen(filename,'w');
 fprintf(fileID,['{"type":"FeatureCollection","features":[{"type":"Feature","properties":{"name":"',char(trackname{k}),'","tags":"Critical_Passage","component":"ocean","author":"Andrew Roberts","object":"transect"},"geometry":{"type":"LineString","coordinates":[']);
 for j=1:length(coords{k}.lats)-1
  fprintf(fileID,['[',num2str(coords{k}.lons(j),'%2.10f'),',',...
                      num2str(coords{k}.lats(j),'%2.10f'),'],']);
 end
 fprintf(fileID,['[',num2str(coords{k}.lons(j),'%2.10f'),',',...
                       num2str(coords{k}.lats(j),'%2.10f'),']']);
 fprintf(fileID,']}}]}');
 fclose(fileID);

 [x,y]=mfwdtran(coords{k}.lats,coords{k}.lons);
 plot(x,y,'r')

end

title('Foxe Basin Through Flow Forced Open for WC14 r03');

ridgepack_fprint('png','E3SM_Foxe_Basin_Passages',1,2)


