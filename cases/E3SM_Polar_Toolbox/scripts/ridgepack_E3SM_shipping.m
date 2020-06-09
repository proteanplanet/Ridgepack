
clear

shiploc='/Users/afroberts/data/SHIPPING';

cd(shiploc)

nwp=true;
%nwp=false;
nsr=true;
%nsr=false;

clf

% plot up
ridgepack_polarm('shipping')

if nwp

 % grab main channels
 nc=ridgepack_clone('Northwest_Passage_Route');

 % break into segments
 k=1;
 kidx(k)=1;
 for i=2:length(nc.latitude.data)
  if ~isnan(nc.latitude.data(i-1)) & isnan(nc.latitude.data(i))
   k=k+1;
   kidx(k)=i;
  end
 end

 % create sectors
 for k=1:length(kidx)-1
  sectornw{k}.lat=nc.latitude.data(kidx(k):kidx(k+1));
  sectornw{k}.lat=sectornw{k}.lat(~isnan(sectornw{k}.lat));
  sectornw{k}.lon=nc.longitude.data(kidx(k):kidx(k+1));
  sectornw{k}.lon=sectornw{k}.lon(~isnan(sectornw{k}.lon));
 end

 % NW passage first ship route
 shiptrack{1}.lat=sectornw{1}.lat;
 shiptrack{1}.lon=sectornw{1}.lon;

 % NW passage second ship route

 % Search for start and and end of segment, split for efficiency
 dists(1:length(sectornw{1}.lat))=100000000000;
 diste(1:length(sectornw{1}.lat))=100000000000;
 for i=1:100:length(sectornw{1}.lat)
  dists(i)=ridgepack_greatcircle(sectornw{2}.lat(1),...
                                 sectornw{2}.lon(1),...
                                 sectornw{1}.lat(i),...
                                 sectornw{1}.lon(i));
  diste(i)=ridgepack_greatcircle(sectornw{2}.lat(end),...
                                 sectornw{2}.lon(end),...
                                 sectornw{1}.lat(i),...
                                 sectornw{1}.lon(i));
 end

 [distmins,sdx]=min(dists);
 [distmine,edx]=min(diste);

 for i=max(sdx-100,1):min(sdx+100,length(sectornw{1}.lat))
  dists(i)=ridgepack_greatcircle(sectornw{2}.lat(1),...
                                 sectornw{2}.lon(1),...
                                 sectornw{1}.lat(i),...
                                 sectornw{1}.lon(i));
 end

 for i=max(edx-100,1):min(edx+100,length(sectornw{1}.lat))
  diste(i)=ridgepack_greatcircle(sectornw{2}.lat(end),...
                                 sectornw{2}.lon(end),...
                                 sectornw{1}.lat(i),...
                                 sectornw{1}.lon(i));
 end

 [distmins,sdx]=min(dists);
 [distmine,edx]=min(diste);

 shiptrack{2}.lat=[sectornw{1}.lat(1:min(sdx,edx))' ...
                   sectornw{2}.lat' ...
                   sectornw{1}.lat(max(sdx,edx):end)'];
 
 shiptrack{2}.lon=[sectornw{1}.lon(1:min(sdx,edx))' ...
                   sectornw{2}.lon' ...
                   sectornw{1}.lon(max(sdx,edx):end)'];
 

 % plot northwest passage
 xcols=colormap(lines(length(sectornw)));
 pattern{1}=':';
 pattern{2}='--';
 for i=1:length(shiptrack)
  [x,y]=mfwdtran(shiptrack{i}.lat,shiptrack{i}.lon);
  plot(x,y,'Color',xcols(i,:),'LineStyle',pattern{i})
 end

 % netcdf northwest
 nc=ridgepack_clone('Northwest_Passage_Route_Variants');

 % break into segments
 k=1;
 kidx(k)=1;
 for i=2:length(nc.latitude.data)
  if ~isnan(nc.latitude.data(i-1)) & isnan(nc.latitude.data(i))
   k=k+1;
   kidx(k)=i;
  end 
 end 

 % create sectors
 for k=1:length(kidx)-1
  sectornwa{k}.lat=nc.latitude.data(kidx(k):kidx(k+1));
  sectornwa{k}.lat=sectornwa{k}.lat(~isnan(sectornwa{k}.lat));
  sectornwa{k}.lon=nc.longitude.data(kidx(k):kidx(k+1));
  sectornwa{k}.lon=sectornwa{k}.lon(~isnan(sectornwa{k}.lon));
 end 

 % create third route
 clear sectornw

 sectornw{1}.lat=shiptrack{1}.lat;
 sectornw{1}.lon=shiptrack{1}.lon;

 sectornw{2}.lat=sectornwa{4}.lat;
 sectornw{2}.lon=sectornwa{4}.lon;

 % Search for start and and end of segment, split for efficiency
 dists(1:length(sectornw{1}.lat))=100000000000;
 diste(1:length(sectornw{1}.lat))=100000000000;
 for i=1:100:length(sectornw{1}.lat)
  dists(i)=ridgepack_greatcircle(sectornw{2}.lat(1),...
                                 sectornw{2}.lon(1),...
                                 sectornw{1}.lat(i),...
                                 sectornw{1}.lon(i));
  diste(i)=ridgepack_greatcircle(sectornw{2}.lat(end),...
                                 sectornw{2}.lon(end),...
                                 sectornw{1}.lat(i),...
                                 sectornw{1}.lon(i));
 end

 [distmins,sdx]=min(dists);
 [distmine,edx]=min(diste);

 for i=max(sdx-100,1):min(sdx+100,length(sectornw{1}.lat))
  dists(i)=ridgepack_greatcircle(sectornw{2}.lat(1),...
                                 sectornw{2}.lon(1),...
                                 sectornw{1}.lat(i),...
                                 sectornw{1}.lon(i));
 end

 for i=max(edx-100,1):min(edx+100,length(sectornw{1}.lat))
  diste(i)=ridgepack_greatcircle(sectornw{2}.lat(end),...
                                 sectornw{2}.lon(end),...
                                 sectornw{1}.lat(i),...
                                 sectornw{1}.lon(i));
 end

 [distmins,sdx]=min(dists);
 [distmine,edx]=min(diste);

 shiptrack{3}.lat=[sectornw{1}.lat(1:min(sdx,edx))' ...
                   sectornw{2}.lat' ...
                   sectornw{1}.lat(max(sdx,edx):end)'];

 shiptrack{3}.lon=[sectornw{1}.lon(1:min(sdx,edx))' ...
                   sectornw{2}.lon' ...
                   sectornw{1}.lon(max(sdx,edx):end)'];

 [x,y]=mfwdtran(shiptrack{3}.lat,shiptrack{3}.lon);
 plot(x,y,'Color','m','LineStyle','-.')

 % create fourth route
 clear sectornw

 sectornw{1}.lat=shiptrack{2}.lat;
 sectornw{1}.lon=shiptrack{2}.lon;

 sectornw{2}.lat=sectornwa{2}.lat';
 sectornw{2}.lon=sectornwa{2}.lon';

 sectornw{3}.lat=shiptrack{1}.lat';
 sectornw{3}.lon=shiptrack{1}.lon';

 % Search for start and and end of segment, split for efficiency
 dists(1:length(sectornw{1}.lat))=100000000000;
 diste(1:length(sectornw{1}.lat))=100000000000;
 for i=1:100:length(sectornw{1}.lat)
  dists(i)=ridgepack_greatcircle(sectornw{2}.lat(1),...
                                 sectornw{2}.lon(1),...
                                 sectornw{1}.lat(i),...
                                 sectornw{1}.lon(i));
 end

 [distmins,sdx]=min(dists);

 for i=max(sdx-100,1):min(sdx+100,length(sectornw{1}.lat))
  dists(i)=ridgepack_greatcircle(sectornw{2}.lat(1),...
                                 sectornw{2}.lon(1),...
                                 sectornw{1}.lat(i),...
                                 sectornw{1}.lon(i));
 end

 [distmins,sdx]=min(dists);

 % Search for start and and end of segment, split for efficiency
 dists(1:length(sectornw{1}.lat))=100000000000;
 diste(1:length(sectornw{1}.lat))=100000000000;
 for i=1:100:length(sectornw{3}.lat)
  diste(i)=ridgepack_greatcircle(sectornw{2}.lat(end),...
                                 sectornw{2}.lon(end),...
                                 sectornw{3}.lat(i),...
                                 sectornw{3}.lon(i));
 end

 [distmine,edx]=min(diste);

 for i=max(edx-100,1):min(edx+100,length(sectornw{3}.lat))
  diste(i)=ridgepack_greatcircle(sectornw{2}.lat(end),...
                                 sectornw{2}.lon(end),...
                                 sectornw{3}.lat(i),...
                                 sectornw{3}.lon(i));
 end

 [distmine,edx]=min(diste);

 shiptrack{4}.lat=[sectornw{1}.lat(1:sdx) ...
                   sectornw{2}.lat ...
                   sectornw{3}.lat(edx:end)];

 shiptrack{4}.lon=[sectornw{1}.lon(1:sdx) ...
                   sectornw{2}.lon ...
                   sectornw{3}.lon(edx:end)];

 [x,y]=mfwdtran(shiptrack{4}.lat,shiptrack{4}.lon);
 plot(x,y,'Color','c','LineStyle','-')

end

if nsr

 % netcdf northern
 nc=ridgepack_clone('Northern_Sea_Route');

 % break into segments
 k=1; 
 kidx(k)=1;
 for i=2:length(nc.latitude.data)
  if ~isnan(nc.latitude.data(i-1)) & isnan(nc.latitude.data(i))
   k=k+1;
   kidx(k)=i;
  end
 end

 % create sectors
 for k=1:length(kidx)-1
  sectorns{k}.lat=nc.latitude.data(kidx(k):kidx(k+1));
  sectorns{k}.lat=sectorns{k}.lat(~isnan(sectorns{k}.lat));
  sectorns{k}.lat=sectorns{k}.lat(end:-1:1);
  sectorns{k}.lon=nc.longitude.data(kidx(k):kidx(k+1));
  sectorns{k}.lon=sectorns{k}.lon(~isnan(sectorns{k}.lon));
  sectorns{k}.lon=sectorns{k}.lon(end:-1:1);
 end

 % now join these together to a single track, one segment ata time
 for j=1:length(sectorns)-2
  dists(1:length(sectorns{j}.lat))=100000000000;
  diste(1:length(sectorns{j}.lat))=100000000000;
  for i=1:100:length(sectorns{j}.lat)
   dists(i)=ridgepack_greatcircle(sectorns{j+1}.lat(1),...
                                  sectorns{j+1}.lon(1),...
                                  sectorns{j}.lat(i),...
                                  sectorns{j}.lon(i));
  end

  [distmins,sdx]=min(dists);

  for i=max(sdx-100,1):min(sdx+100,length(sectorns{j}.lat))
   dists(i)=ridgepack_greatcircle(sectorns{j+1}.lat(1),...
                                  sectorns{j+1}.lon(1),...
                                  sectorns{j}.lat(i),...
                                  sectorns{j}.lon(i));
  end

  [distmins,sidx(j)]=min(dists);

 end

 shiptrack{5}.lat=[sectorns{1}.lat(1:sidx(1))' ...
                   sectorns{2}.lat(1:sidx(2))' ...
                   sectorns{3}.lat(1:end)'];
 
 shiptrack{5}.lon=[sectorns{1}.lon(1:sidx(1))' ...
                   sectorns{2}.lon(1:sidx(2))' ...
                   sectorns{3}.lon(1:end)'];
 
 % plot northern sea route
 [x,y]=mfwdtran(shiptrack{3}.lat,shiptrack{3}.lon);
 plot(x,y,'Color','b','LineStyle','--')

 % netcdf northern
 nc=ridgepack_clone('Northern_Sea_Route_Variants');

 % break into segments
 k=1;
 kidx(k)=1;
 for i=2:length(nc.latitude.data)

  if ~isnan(nc.latitude.data(i-1)) & isnan(nc.latitude.data(i))
   k=k+1;
   kidx(k)=i;
  end

 end

 % create sector
 for k=1:length(kidx)-1
  sectornsa{k}.lat=nc.latitude.data(kidx(k):kidx(k+1));
  sectornsa{k}.lat=sectornsa{k}.lat(~isnan(sectornsa{k}.lat));
  sectornsa{k}.lon=nc.longitude.data(kidx(k):kidx(k+1));
  sectornsa{k}.lon=sectornsa{k}.lon(~isnan(sectornsa{k}.lon));
 end

 clear sectorns

 sectorns{1}.lat=shiptrack{3}.lat;
 sectorns{1}.lon=shiptrack{3}.lon;

 sectorns{2}.lat=[sectornsa{5}.lat' sectornsa{6}.lat'];
 sectorns{2}.lon=[sectornsa{5}.lon' sectornsa{6}.lon'];

 % Search for start and and end of segment, split for efficiency
 dists(1:length(sectorns{1}.lat))=100000000000;
 diste(1:length(sectorns{1}.lat))=100000000000;
 for i=1:100:length(sectorns{1}.lat)
  dists(i)=ridgepack_greatcircle(sectorns{2}.lat(1),...
                                 sectorns{2}.lon(1),...
                                 sectorns{1}.lat(i),...
                                 sectorns{1}.lon(i));
  diste(i)=ridgepack_greatcircle(sectorns{2}.lat(end),...
                                 sectorns{2}.lon(end),...
                                 sectorns{1}.lat(i),...
                                 sectorns{1}.lon(i));
 end

 [distmins,sdx]=min(dists);
 [distmine,edx]=min(diste);

 for i=max(sdx-100,1):min(sdx+100,length(sectorns{1}.lat))
  dists(i)=ridgepack_greatcircle(sectorns{2}.lat(1),...
                                 sectorns{2}.lon(1),...
                                 sectorns{1}.lat(i),...
                                 sectorns{1}.lon(i));
 end

 for i=max(edx-100,1):min(edx+100,length(sectorns{1}.lat))
  diste(i)=ridgepack_greatcircle(sectorns{2}.lat(end),...
                                 sectorns{2}.lon(end),...
                                 sectorns{1}.lat(i),...
                                 sectorns{1}.lon(i));
 end

 [distmins,sdx]=min(dists);
 [distmine,edx]=min(diste);

 shiptrack{6}.lat=[sectorns{1}.lat(1:min(sdx,edx)) ...
                   sectorns{2}.lat(end:-1:1) ...
                   sectorns{1}.lat(max(sdx,edx):end)];

 shiptrack{6}.lon=[sectorns{1}.lon(1:min(sdx,edx)) ...
                   sectorns{2}.lon(end:-1:1) ...
                   sectorns{1}.lon(max(sdx,edx):end)];

 % plot northern sea route
 [x,y]=mfwdtran(shiptrack{4}.lat,shiptrack{4}.lon);
 plot(x,y,'Color','r','LineStyle',':')

end

