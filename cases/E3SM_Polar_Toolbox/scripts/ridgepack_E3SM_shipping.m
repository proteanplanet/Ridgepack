
clear

shiploc='/Users/afroberts/data/SHIPPING';

nwp=true;
nsr=true;

%nwp=false;
%nsr=false;

%discretize=true;
discretize=false;

clf

cd(shiploc)

% plot up
ridgepack_polarm('shipping')

% plot northwest passage
xcols=colormap(lines(6));

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
 clear dists diste
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
 
 [x,y]=mfwdtran(shiptrack{1}.lat,shiptrack{1}.lon);
 plot(x,y,'Color',xcols(1,:),'LineStyle','-')

 [x,y]=mfwdtran(shiptrack{2}.lat,shiptrack{2}.lon);
 plot(x,y,'Color',xcols(2,:),'LineStyle','--')

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
 clear dists diste
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
 clear dists diste
 dists(1:length(sectornw{1}.lat))=100000000000;
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

 % rearrange NWP tracks to move from Pacific to Atlantic
 for i=1:4
  lats=shiptrack{i}.lat(end:-1:1);
  lons=shiptrack{i}.lon(end:-1:1);
  shiptrack{i}.lat=lats;
  shiptrack{i}.lon=lons;
 end

end

drawnow

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

 % now join these together to a single track, one segment at a time
 for j=1:length(sectorns)-2
  clear dists
  dists(1:length(sectorns{j}.lat))=100000000000;
  for i=1:100:length(sectorns{j}.lat)
   dists(i)=ridgepack_greatcircle(sectorns{j+1}.lat(1),...
                                  sectorns{j+1}.lon(1),...
                                  sectorns{j}.lat(i),...
                                  sectorns{j}.lon(i));
  end
  [distmins,sdx]=min(dists);

  clear dists
  dists(1:min(sdx+100,length(sectorns{j}.lat)))=100000000000;
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
                   sectorns{3}.lat(30:end)'];
 
 shiptrack{5}.lon=[sectorns{1}.lon(1:sidx(1))' ...
                   sectorns{2}.lon(1:sidx(2))' ...
                   sectorns{3}.lon(30:end)'];

 % join up with NWP 
 disti(1:length(shiptrack{1}.lat))=100000000000;
 distj(1:length(shiptrack{5}.lat))=100000000000;

 for i=1:100:length(shiptrack{1}.lat)
 for j=1:10:length(shiptrack{5}.lat)
  distj(j)=ridgepack_greatcircle(shiptrack{1}.lat(i),...
                                 shiptrack{1}.lon(i),...
                                 shiptrack{5}.lat(j),...
                                 shiptrack{5}.lon(j));
 end
 [disti(i),jdx(i)]=min(distj);
 distj(1:length(shiptrack{5}.lat))=100000000000;
 end
 II=find(disti==min(disti));

 disti(1:min(II+100,length(shiptrack{1}.lat)))=100000000000;
 distj(1:min(jdx(II)+10,length(shiptrack{5}.lat)))=100000000000;

 for i=max(1,II-100):min(II+100,length(shiptrack{1}.lat))
 for j=max(1,jdx(II)-10):min(jdx(II)+10,length(shiptrack{5}.lat))
  distj(j)=ridgepack_greatcircle(shiptrack{1}.lat(i),...
                                 shiptrack{1}.lon(i),...
                                 shiptrack{5}.lat(j),...
                                 shiptrack{5}.lon(j));
 end
 [disti(i),jdx(i)]=min(distj);
 distj(1:min(jdx(II)+10,length(shiptrack{5}.lat)))=100000000000;
 end
 II=find(disti==min(disti))
 JJ=jdx(II)

 shiplat=[shiptrack{1}.lat(end:-1:II)' shiptrack{5}.lat(JJ:end)];
 shiplon=[shiptrack{1}.lon(end:-1:II)' shiptrack{5}.lon(JJ:end)];

 shiptrack{5}.lat=shiplat;
 shiptrack{5}.lon=shiplon;
 
 % plot northern sea route
 [x,y]=mfwdtran(shiptrack{5}.lat,shiptrack{5}.lon);
 plot(x,y,'Color',xcols(5,:),'LineStyle','-')

 drawnow

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

 % generate sixth ship track
 clear sectorns

 sectorns{1}.lat=shiptrack{5}.lat;
 sectorns{1}.lon=shiptrack{5}.lon;

 sectorns{2}.lat=[sectornsa{5}.lat' sectornsa{6}.lat'];
 sectorns{2}.lon=[sectornsa{5}.lon' sectornsa{6}.lon'];

 % Search for start and and end of segment, split for efficiency
 clear dists diste
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

 shiptrack{6}.lat=[sectorns{1}.lat(1:edx) ...
                   sectorns{2}.lat(end:-1:1) ...
                   sectorns{1}.lat(sdx:end)];

 shiptrack{6}.lon=[sectorns{1}.lon(1:edx) ...
                   sectorns{2}.lon(end:-1:1) ...
                   sectorns{1}.lon(sdx:end)];

 % plot northern sea route
 [x,y]=mfwdtran(shiptrack{6}.lat,shiptrack{6}.lon);
 plot(x,y,'Color',xcols(6,:),'LineStyle','-')

end

if nwp | nsr

 save('total_tracks','shiptrack')

else

 load('total_tracks')

 % plot tracks
 ridgepack_polarm('shipping') 
 xcols=colormap(lines(length(shiptrack)));
 for i=1:length(shiptrack)
  [x,y]=mfwdtran(shiptrack{i}.lat,shiptrack{i}.lon);
  plot(x,y,'Color',xcols(i,:),'LineStyle','-')
 end
 ridgepack_fprint('png','test',1,2)

end

% turn the data into 1 NM - seperated points
if discretize

 %for i=1:1
 for i=1:length(shiptrack)
  k=1;
  lats(k)=shiptrack{i}.lat(1);
  lons(k)=shiptrack{i}.lon(1);
  for j=2:length(shiptrack{i}.lat)
   theta=pi/180;
   while theta*180/pi>1/60
    k=k+1;
    [x,y,z,phi,theta]=ridgepack_satfwd(shiptrack{i}.lat(j),...
                                       shiptrack{i}.lon(j),...
                                       lats(k-1),lons(k-1),90,1,0);
    [lats(k),lons(k)]=ridgepack_satinv(phi,(pi/180)*(1/60),...
                                       lats(k-1),lons(k-1));
   end
  end

  % check data to make sure there is 1 NM spacing
  theta=zeros([1 length(lats)-1]);
  for k=2:length(lats)
   [x,y,z,phi,theta(k-1)]=ridgepack_satfwd(lats(k),lons(k),...
                                           lats(k-1),lons(k-1),...
                                           90,1,0);
  end
  if any(round(theta(:)*60*180/pi,6)~=1)
   error('not NM compliant')
  end
  [xp,yp]=mfwdtran(lats,lons);
  plot(xp,yp,'Color',xcols(i,:))
  drawnow
  clear lats lons
 end

end

return

clf


%annotation(figure1,'line',[0.610714285714286 0.271428571428571],...
%    [0.828571428571429 0.457142857142857]);

ridgepack_polarm('shipping')

xp=[0.610714285714286 0.271428571428571];
yp=[0.828571428571429 0.457142857142857];
AX=axis(gca);
Xrange=AX(2)-AX(1);
Yrange=AX(4)-AX(3);       
PO=get(gca,'Position');

x(1:2)=AX(1)+Xrange*(xp(1:2)-PO(1))/PO(3);
y(1:2)=AX(3)+Yrange*(yp(1:2)-PO(2))/PO(4);

[lats lons]=minvtran(x,y);

[dist,angl,phi,tracklat,tracklon,tracklen]=...
          ridgepack_greatcircle(lats(1),lons(1),lats(2),lons(2));

[x,y]=mfwdtran(tracklat,tracklon);

plot(x,y,'r')


% Create line
%annotation(figure1,'line',[0.275 0.733928571428571],...
%    [0.327571428571429 0.721428571428572]);

% Create line
%annotation(figure1,'line',[0.289285714285714 0.655357142857143],...
%    [0.402380952380952 0.761904761904762]);



