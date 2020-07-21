
clear

shiploc='/Users/afroberts/data/SHIPPING';

%generate=true;
generate=false;

%discretize=true;
discretize=false;

%newnw=true;
newnw=false;

%e3sm=true;
e3sm=false;

clf

cd(shiploc)

% plot up
bering=false;
%bering=true;
if bering
 ridgepack_conicm('bering')
else
 ridgepack_polarm('shipping')
end

% plot northwest passage
xcols=colormap(lines(20));

if generate

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

 % make conglomerate track of track 6 and track 8
 xs=length(shiptrack{1}.lat)-42000;
 xe=length(shiptrack{1}.lat)-30000;
 [x,y,z,phi,theta]=ridgepack_satfwd(shiptrack{1}.lat(xe),...
                                    shiptrack{1}.lon(xe),...
				    shiptrack{1}.lat(xs),...
                                    shiptrack{1}.lon(xs),...
                                    90,1,0);

 clear lat lon

 [lat(1),lon(1)]=ridgepack_satinv(0.78*phi,0.02*theta,...
 			          shiptrack{1}.lat(xs),...
                                  shiptrack{1}.lon(xs));
 maxbend=50
 for i=2:maxbend/2
  [lat(i),lon(i)]=ridgepack_satinv((0.78+0.45*i/maxbend)*phi,0.02*theta,...
                                  lat(i-1),lon(i-1));
 end

 ml=length(lat);
 xe=length(shiptrack{1}.lat)-27700;

 % pad lat and lon so that final length of shiptrack is
 % same as original. This is needed because fix to Bering Strait
 % was a late hange
 for i=ml+1:xe-xs+1
  lat(i)=lat(ml);
  lon(i)=lon(ml);
 end

 shiptrack{1}.lat=[shiptrack{1}.lat(1:xs)' lat ...
                   shiptrack{1}.lat(xe:end)'];
 shiptrack{1}.lon=[shiptrack{1}.lon(1:xs)' lon ...
                   shiptrack{1}.lon(xe:end)'];

 sectornw{1}.lat=shiptrack{1}.lat;
 sectornw{1}.lon=shiptrack{1}.lon;

 clear lat lon

 [x,y]=mfwdtran(shiptrack{1}.lat,shiptrack{1}.lon);
 plot(x,y,'Color',xcols(1,:),'LineStyle','-')
 shiptrack{1}.name='Northwest Passage: M''Clintock Channel Variant';

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

 max(sdx,edx)

 shiptrack{2}.lat=[sectornw{1}.lat(1:min(sdx,edx)) ...
                   sectornw{2}.lat' NaN ...
                   sectornw{1}.lat(219000:end)];
 
 shiptrack{2}.lon=[sectornw{1}.lon(1:min(sdx,edx)) ...
                   sectornw{2}.lon' NaN ...
                   sectornw{1}.lon(219000:end)];
 
 [x,y]=mfwdtran(shiptrack{2}.lat,shiptrack{2}.lon);
 %plot(x,y,'Color',xcols(2,:),'LineStyle','-')
 shiptrack{2}.name='Northwest Passage: M''Clure Strait Variant';

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

 shiptrack{3}.lat=[sectornw{1}.lat(1:min(sdx,edx)) ...
                   sectornw{2}.lat' ...
                   sectornw{1}.lat(max(sdx,edx):end)];

 shiptrack{3}.lon=[sectornw{1}.lon(1:min(sdx,edx)) ...
                   sectornw{2}.lon' ...
                   sectornw{1}.lon(max(sdx,edx):end)];

 [x,y]=mfwdtran(shiptrack{3}.lat,shiptrack{3}.lon);
 plot(x,y,'Color','m','LineStyle','-')
 shiptrack{3}.name='Northwest Passage: Peel Sound Variant';

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

 edx=213404;

 shiptrack{4}.lat=[sectornw{1}.lat(1:sdx) ...
                   sectornw{2}.lat ...
                   sectornw{3}.lat(edx:end)'];

 shiptrack{4}.lon=[sectornw{1}.lon(1:sdx) ...
                   sectornw{2}.lon ...
                   sectornw{3}.lon(edx:end)'];

 [x,y]=mfwdtran(shiptrack{4}.lat,shiptrack{4}.lon);
 plot(x,y,'Color','b','LineStyle','-')
 shiptrack{4}.name='Northwest Passage: Prince of Wales Strait Variant';

 % rearrange NWP tracks to move from Pacific to Atlantic
 for i=1:4
  lats=shiptrack{i}.lat(end:-1:1);
  lons=shiptrack{i}.lon(end:-1:1);
  shiptrack{i}.lat=lats;
  shiptrack{i}.lon=lons;
 end

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
  if j==1
   pointss=1;
  elseif j==2 
   pointss=15000;
  end
  clear dists
  dists(1:length(sectorns{j}.lat))=100000000000;
  for i=1:1000:length(sectorns{j}.lat)
   dists(i)=ridgepack_greatcircle(sectorns{j+1}.lat(pointss),...
                                  sectorns{j+1}.lon(pointss),...
                                  sectorns{j}.lat(i),...
                                  sectorns{j}.lon(i));
  end
  [distmins,sdx]=min(dists);
 
  clear dists
  dists(1:min(sdx+1000,length(sectorns{j}.lat)))=100000000000;
  for i=max(sdx-1000,1):min(sdx+1000,length(sectorns{j}.lat))
   dists(i)=ridgepack_greatcircle(sectorns{j+1}.lat(pointss),...
                                  sectorns{j+1}.lon(pointss),...
                                  sectorns{j}.lat(i),...
                                  sectorns{j}.lon(i));
  end
  [distmins,sidx(j)]=min(dists);

 end

 shiptrack{5}.lat=[sectorns{1}.lat(1:sidx(1))' ...
                   sectorns{2}.lat(1:sidx(2))' ...
                   sectorns{3}.lat(16000:end)'];
 
 shiptrack{5}.lon=[sectorns{1}.lon(1:sidx(1))' ...
                   sectorns{2}.lon(1:sidx(2))' ...
                   sectorns{3}.lon(16000:end)'];

 shiplat=[shiptrack{1}.lat(1:45000) shiptrack{5}.lat(50:end)];
 shiplon=[shiptrack{1}.lon(1:45000) shiptrack{5}.lon(50:end)];

 shiptrack{5}.lat=shiplat;
 shiptrack{5}.lon=shiplon;
 
 % plot northern sea route
 [x,y]=mfwdtran(shiptrack{5}.lat,shiptrack{5}.lon);
 plot(x,y,'Color',xcols(5,:),'LineStyle','-')
 shiptrack{5}.name='Northern Sea Route: Wrangel Island Detour';

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
 shiptrack{6}.name='Northern Sea Route: Proliv Longa Variant';

 % now add in Great Circle Route
 nc=ridgepack_clone('Great_Circle_Route');
 shiptrack{7}.lat=[shiptrack{1}.lat(1:52000) ...
                   nc.latitude.data(700:-1:1)' ...
                   shiptrack{6}.lat(end)];
 shiptrack{7}.lon=[shiptrack{1}.lon(1:52000) ...
                   nc.longitude.data(700:-1:1)' ...
                   shiptrack{6}.lon(end)];

 [x,y]=mfwdtran(shiptrack{7}.lat,shiptrack{7}.lon);
 plot(x,y,'Color',xcols(7,:),'LineStyle','-')
 shiptrack{7}.name='Arctic Ocean Great Circle Route';
 
 drawnow

 ridgepack_fprint('png','raw_tracks',1,2)

 save('raw_tracks','shiptrack')

elseif discretize

 load('raw_tracks')

 % turn the data into 1 NM - seperated points
 ridgepack_polarm('shipping') 
 
 for i=1:length(shiptrack)
   k=1;
   lats(k)=shiptrack{i}.lat(1);
   lons(k)=shiptrack{i}.lon(1);
   for j=2:length(shiptrack{i}.lat)
    theta=pi/180;
    while theta*180/pi>1/60
     [x,y,z,phi,theta]=ridgepack_satfwd(shiptrack{i}.lat(j),...
                                        shiptrack{i}.lon(j),...
                                        lats(k),lons(k),90,1,0);
     if theta*180/pi>=1/60
      k=k+1;
      [lats(k),lons(k)]=ridgepack_satinv(phi,(pi/180)*(1/60),...
                                        lats(k-1),lons(k-1));
     end

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
    error('not nautical mile compliant')
   end

   % create distance along track in nautical miles
   % and fraction along track
   for k=1:length(lats)
    distalt(k)=k-1;
    fractiondist(k)=distalt(k)/(length(lats)-1);
   end

   % plot the data
   [xp,yp]=mfwdtran(lats,lons);
   plot(xp,yp,'.','Color',xcols(i,:))

   drawnow

   track{i}.lat=lats;
   track{i}.latunits='degrees north';
   track{i}.lon=lons;
   track{i}.lonunits='degrees east';
   track{i}.name=char(shiptrack{i}.name);
   track{i}.distance=distalt;
   track{i}.distanceunits='nautical miles';
   track{i}.fractionaldistance=fractiondist;

   clear lats lons distalt fractiondist

 end

 ridgepack_fprint('png','interim_tracks',1,2)

 save('interim_tracks','track')

% build new tracks based on combinations
elseif newnw

 load('interim_tracks')

 % clean up first made-up track 
 for i=1:length(track{7}.lat) 

  diffdist=track{5}.fractionaldistance-track{7}.fractionaldistance(i);

  idx=find(abs(diffdist)==min(abs(diffdist)));

  [x,y,z,phi,theta]=ridgepack_satfwd(track{5}.lat(idx),...
                                     track{5}.lon(idx),...
				     track{7}.lat(i),...
                                     track{7}.lon(i),...
                                     90,1,0);

  [latsa(i),lonsa(i)]=ridgepack_satinv(phi,0.70*theta,...
	  			       track{7}.lat(i),...
                                       track{7}.lon(i));
  
  [latsb(i),lonsb(i)]=ridgepack_satinv(phi,0.43*theta,...
	  			       track{7}.lat(i),...
                                       track{7}.lon(i));
  
  [latsc(i),lonsc(i)]=ridgepack_satinv(phi,0.16*theta,...
	  			       track{7}.lat(i),...
                                       track{7}.lon(i));
  
 end

 [x,y,z,phi,theta]=ridgepack_satfwd(track{7}.lat(end),...
                                    track{7}.lon(end),...
                                    latsa(1995),...
                                    lonsa(1995),...
                                    90,1,0);

 [lat(1),lon(1)]=ridgepack_satinv(0.82*phi,0.02*theta,...
                                  latsa(2000),lonsa(2000));
 for i=2:26
  [lat(i),lon(i)]=ridgepack_satinv(0.82*phi,0.02*theta,...
                                  lat(i-1),lon(i-1));
 end

 lat(end+1)=latsa(end);
 lon(end+1)=lonsa(end);

 shiptrack{10}.lat=[track{5}.lat(1:900) latsa(800:1200) ...
                    latsa(1300:1400) latsa(1700:2000) lat];
 shiptrack{10}.lon=[track{5}.lon(1:900) lonsa(800:1200) ...
                    lonsa(1300:1400) lonsa(1700:2000) lon];
 shiptrack{10}.name='Arctic Ocean Route: Franz Josef Land Detour';

 clear lat lon

 % clean up second made-up track 
 [x,y,z,phi,theta]=ridgepack_satfwd(track{7}.lat(end),...
                                    track{7}.lon(end),...
                                    latsb(2150),...
                                    lonsb(2150),...
                                    90,1,0);

 [lat(1),lon(1)]=ridgepack_satinv(0.84*phi,0.02*theta,...
                                  latsb(2150),lonsb(2150));
 for i=2:35
  [lat(i),lon(i)]=ridgepack_satinv(0.84*phi,0.02*theta,...
                                  lat(i-1),lon(i-1));
 end

 lat(end+1)=latsb(end);
 lon(end+1)=lonsb(end);

 shiptrack{9}.lat=[track{5}.lat(1:700) latsb(800:1100) ...
                   latsb(1300:1350) latsb(2148:2150) lat];
 shiptrack{9}.lon=[track{5}.lon(1:700) lonsb(800:1100) ...
                   lonsb(1300:1350) lonsb(2148:2150) lon];
 shiptrack{9}.name='Arctic Ocean Route: Nordaustlandet Detour';

 clear lat lon

 % clean up third made-up track 
 [x,y,z,phi,theta]=ridgepack_satfwd(track{7}.lat(end),...
                                    track{7}.lon(end),...
                                    latsc(2250),...
                                    lonsc(2250),...
                                    90,1,0);

 [lat(1),lon(1)]=ridgepack_satinv(0.94*phi,0.02*theta,...
                                  latsc(2250),lonsc(2250));
 for i=2:40
  [lat(i),lon(i)]=ridgepack_satinv(0.94*phi,0.02*theta,...
                                  lat(i-1),lon(i-1));
 end

 lat(end+1)=latsc(end);
 lon(end+1)=lonsc(end);

 shiptrack{8}.lat=[track{5}.lat(1:560) latsc(800:1200) ...
                   latsc(2200:2250) lat];
 shiptrack{8}.lon=[track{5}.lon(1:560) lonsc(800:1200) ...
                   lonsc(2200:2250) lon];
 shiptrack{8}.name='Arctic Ocean Route: Fram Strait Detour 1';

 % clean up track 8 further using google Earth help
 kmlStructT=kml2struct('Transpolar_Variation.kml');

 shiptrack{8}.lat=[shiptrack{8}.lat(1:560)  ...
                   kmlStructT.Lat(1:end-1)' ...
                   shiptrack{8}.lat(end)];
 shiptrack{8}.lon=[shiptrack{8}.lon(1:560)  ...
                   kmlStructT.Lon(1:end-1)' ...
                   shiptrack{8}.lon(end)];

 clear lat lon

 % make conglomerate track of track 6 and track 8
 [x,y,z,phi,theta]=ridgepack_satfwd(track{6}.lat(end),...
                                    track{6}.lon(end),...
                                    track{6}.lat(2350),...
                                    track{6}.lon(2350),...
                                    90,1,0);

 [lat(1),lon(1)]=ridgepack_satinv(1.15*phi,0.01*theta,...
                                  track{6}.lat(2350),...
                                  track{6}.lon(2350));
 for i=2:20
  [lat(i),lon(i)]=ridgepack_satinv(1.15*phi,0.01*theta,...
                                  lat(i-1),lon(i-1));
 end

 for i=21:24
  [lat(i),lon(i)]=ridgepack_satinv(1.20*phi,0.01*theta,...
                                  lat(i-1),lon(i-1));
 end

 for i=25:28
  [lat(i),lon(i)]=ridgepack_satinv(1.22*phi,0.01*theta,...
                                  lat(i-1),lon(i-1));
 end

 shiptrack{11}.lat=[track{6}.lat(1:2350) lat ...
                    shiptrack{10}.lat(1720:end)];
 shiptrack{11}.lon=[track{6}.lon(1:2350) lon ...
                    shiptrack{10}.lon(1720:end)];
 shiptrack{11}.name='Northern Sea Route: Kara Shortcut with Proliv Longa Passage';

 clear lat lon

 % make conglomerate track of track 5 and track 8
 [x,y,z,phi,theta]=ridgepack_satfwd(track{5}.lat(end),...
                                    track{5}.lon(end),...
                                    track{5}.lat(2350),...
                                    track{5}.lon(2350),...
                                    90,1,0);

 [lat(1),lon(1)]=ridgepack_satinv(1.15*phi,0.01*theta,...
                                  track{5}.lat(2350),...
                                  track{5}.lon(2350));
 for i=2:20
  [lat(i),lon(i)]=ridgepack_satinv(1.15*phi,0.01*theta,...
                                  lat(i-1),lon(i-1));
 end

 for i=21:24
  [lat(i),lon(i)]=ridgepack_satinv(1.20*phi,0.01*theta,...
                                  lat(i-1),lon(i-1));
 end

 for i=25:28
  [lat(i),lon(i)]=ridgepack_satinv(1.22*phi,0.01*theta,...
                                  lat(i-1),lon(i-1));
 end

 shiptrack{12}.lat=[track{5}.lat(1:2350) lat ...
                    shiptrack{10}.lat(1720:end)];
 shiptrack{12}.lon=[track{5}.lon(1:2350) lon ...
                    shiptrack{10}.lon(1720:end)];
 shiptrack{12}.name='Northern Sea Route: Kara Shortcut and Wrangel Island Detour';

 clear lat lon

 % make conglomerate of track 5 and 10
 kmlStructA=kml2struct('ConnectorA.kml');

 shiptrack{13}.lat=[track{5}.lat(1:1640) ...
                    kmlStructA.Lat(end:-1:1)' ...
                    shiptrack{10}.lat(1460:end)];
 shiptrack{13}.lon=[track{5}.lon(1:1640) ...
                    kmlStructA.Lon(end:-1:1)' ...
                    shiptrack{10}.lon(1460:end)];
 shiptrack{13}.name='Northern Sea Route: Siberian Coastal Route with Franz Josef Land Shortcut';

 clear lat lon

 % make conglomerate of track 5 and 9 
 kmlStructB=kml2struct('ConnectorB.kml');

 shiptrack{14}.lat=[track{5}.lat(1:1640) ...
                    kmlStructA.Lat(end:-1:1)' ...
                    kmlStructB.Lat(end:-1:1)' ...
                    shiptrack{9}.lat(1061:end)];
 shiptrack{14}.lon=[track{5}.lon(1:1640) ...
                    kmlStructA.Lon(end:-1:1)' ...
                    kmlStructB.Lon(end:-1:1)' ...
                    shiptrack{9}.lon(1061:end)];
 shiptrack{14}.name='Northern Sea Route: Siberian Coastal Route with Svalbard Shortcut';

 clear lat lon

 % make conglomerate of track 6 and 10
 shiptrack{15}.lat=[track{6}.lat(1:1040) ...
                    shiptrack{13}.lat(1080:end)];
 shiptrack{15}.lon=[track{6}.lon(1:1040) ...
                    shiptrack{13}.lon(1080:end)];
 shiptrack{15}.name='Northern Sea Route: Siberian Coastal Route inside Wrangel Island with Svalbard Shortcut';
                    
 clear lat lon

 % make conglomerate of track 6 and 9 
 shiptrack{16}.lat=[track{6}.lat(1:1040) ...
                    shiptrack{14}.lat(1080:end)];
 shiptrack{16}.lon=[track{6}.lon(1:1040) ...
                    shiptrack{14}.lon(1080:end)];
 shiptrack{16}.name='Northern Sea Route: Siberian Coastal Route inside Wrangel Island with Svalbard Shortcut';
                    
 clear lat lon

 % make conglomerate of track 9 and 10
 shiptrack{17}.lat=[shiptrack{10}.lat(1:1460) ...
                    kmlStructB.Lat(end:-1:1)' ...
                    shiptrack{9}.lat(1061:end)];
 shiptrack{17}.lon=[shiptrack{10}.lon(1:1460) ...
                    kmlStructB.Lon(end:-1:1)' ...
                    shiptrack{9}.lon(1061:end)];
 shiptrack{17}.name='Northern Sea Route: Siberian Coastal Route with Svalbard Shortcut';

 clear lat lon

 % add in further variant of track 9 to go through Frame Strait
 kmlStructC=kml2struct('ConnectorC.kml');

 shiptrack{18}.lat=[shiptrack{9}.lat(1:1050)  ...
                   kmlStructC.Lat(end-1:-1:2)' ...
                   shiptrack{9}.lat(end)];
 shiptrack{18}.lon=[shiptrack{9}.lon(1:1050)  ...
                   kmlStructC.Lon(end-1:-1:2)' ...
                   shiptrack{9}.lon(end)];
 shiptrack{18}.name='Arctic Ocean Route: Fram Strait Detour 2';

 clear lat lon

 % add in further variant of track 10 to go through Frame Strait
 kmlStructD=kml2struct('ConnectorD.kml');

 shiptrack{19}.lat=[shiptrack{10}.lat(1:1350)  ...
                   kmlStructD.Lat(end-1:-1:2)' ...
                   shiptrack{18}.lat(1058:end)];
 shiptrack{19}.lon=[shiptrack{10}.lon(1:1350)  ...
                   kmlStructD.Lon(end-1:-1:2)' ...
                   shiptrack{18}.lon(1058:end)];
 shiptrack{19}.name='Arctic Ocean Route: Fram Strait Detour 3';

 clear lat lon


 % plot up the tracks and seperate by 1nm
 ridgepack_polarm('shipping') 

 xdf=[8:19];
 %xdf=10;
 
 for i=xdf

   k=1;
   lats(k)=shiptrack{i}.lat(1);
   lons(k)=shiptrack{i}.lon(1);
   for j=2:length(shiptrack{i}.lat)
    theta=pi/180;
    while theta*180/pi>1/60
     [x,y,z,phi,theta]=ridgepack_satfwd(shiptrack{i}.lat(j),...
                                        shiptrack{i}.lon(j),...
                                        lats(k),lons(k),90,1,0);
     if theta*180/pi>=1/60
      k=k+1;
      [lats(k),lons(k)]=ridgepack_satinv(phi,(pi/180)*(1/60),...
                                        lats(k-1),lons(k-1));
     end

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
    error('not nautical mile compliant')
   end

   % create distance along track in nautical miles
   % and fraction along track
   for k=1:length(lats)
    distalt(k)=k-1;
    fractiondist(k)=distalt(k)/(length(lats)-1);
   end

   % plot the data
   [xp,yp]=mfwdtran(lats,lons);
   plot(xp,yp,'-','Color',xcols(i,:))

   drawnow

   track{i}.lat=lats;
   track{i}.latunits='degrees north';
   track{i}.lon=lons;
   track{i}.lonunits='degrees east';
   track{i}.name=char(shiptrack{i}.name);
   track{i}.distance=distalt;
   track{i}.distanceunits='nautical miles';
   track{i}.fractionaldistance=fractiondist;

   clear lats lons distalt fractiondist

 end

 for i=1:7
  [x,y]=mfwdtran(track{i}.lat,track{i}.lon);
  plot(x,y,'Color',xcols(i,:),'LineStyle','-')
 end

 drawnow

 ridgepack_fprint('png','final_tracks',1,2)

 save('final_tracks','track')


elseif e3sm

  % generate E3SM mesh KML file for mesh generation
  load('final_tracks')

  clf

  ridgepack_polarm('shipping') 

  clear coords trackname

  k=0;

  for i=[1:4 6]

   if i==1

    k=k+1;
    trackname{k}='Bering_Strait_Ship_Track';
    coords{k}.lats=track{i}.lat(150:600);
    coords{k}.lons=track{i}.lon(150:600);

    k=k+1;
    trackname{k}='Northwest_Passage_McClintock_Channel_Variant';
    coords{k}.lats=track{i}.lat(1400:3000);
    coords{k}.lons=track{i}.lon(1400:3000);

   elseif i==2

    k=k+1;
    trackname{k}='Northwest_Passage_McClure_Strait_Variant';
    coords{k}.lats=track{i}.lat(1510:2600);
    coords{k}.lons=track{i}.lon(1510:2600);

   elseif i==3

    k=k+1;
    trackname{k}='Northwest_Passage_Peel_Sound_Variant';
    coords{k}.lats=track{i}.lat(2200:2600);
    coords{k}.lons=track{i}.lon(2200:2600);

   elseif i==4

    k=k+1;
    trackname{k}='Northwest_Passage_Prince_of_Wales_Strait_Variant';
    coords{k}.lats=track{i}.lat(1550:1880);
    coords{k}.lons=track{i}.lon(1550:1880);

   elseif i==6

    k=k+1;
    trackname{k}='Northern_Sea_Route_Proliv_Longa';
    coords{k}.lats=track{i}.lat(750:950);
    coords{k}.lons=track{i}.lon(750:950);
    idx=find(coords{k}.lons>0);
    coords{k}.lats=coords{k}.lats(idx);
    coords{k}.lons=coords{k}.lons(idx);

    k=k+1;
    trackname{k}='Northern_Sea_Route_New_Siberian_Islands';
    coords{k}.lats=track{i}.lat(1450:1750);
    coords{k}.lons=track{i}.lon(1450:1750);

    k=k+1;
    trackname{k}='Northern_Sea_Route_Severnaya_Zemlya';
    coords{k}.lats=track{i}.lat(2050:2350);
    coords{k}.lons=track{i}.lon(2050:2350);

    k=k+1;
    trackname{k}='Northern_Sea_Route_Novaya_Zemlya_Passage';
    coords{k}.lats=track{i}.lat(2910:3290);
    coords{k}.lons=track{i}.lon(2910:3290);

   end

  end


  for k=1:length(trackname)

   % write out GeoJSON
   filename=['E3SM_Shipping_Channel_',num2str(k,'%2.2i'),'.geojson'];
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
   plot(x,y,'b')

  end

  title('Arctic Shipping Passages Forced Open for WC14 r03');

  ridgepack_fprint('png','E3SM_Forced_Arctic_Passages',1,2)

else

  load('final_tracks')

  % plot tracks
  ridgepack_polarm('shipping') 
  xcols=colormap(lines(length(track)));

  nc.attributes.title='InteRFACE Arctic Ship Track Variants';

  for i=1:length(track)

   [x,y]=mfwdtran(track{i}.lat,track{i}.lon);
   plot(x,y,'Color',xcols(i,:),'LineStyle','-')

   name=['track',num2str(i,'%2.2i')];

   nc.([name,'_length']).data=[1:length(track{i}.lat)];
   nc.([name,'_length']).dimension={[name,'_length']};
   nc.([name,'_length']).long_name=[name,' index'];

   nc.([name,'_latitude']).data=track{i}.lat;
   nc.([name,'_latitude']).units='degrees_north';
   nc.([name,'_latitude']).dimension={[name,'_length']};
   nc.([name,'_latitude']).long_name=['latitude, ',...
                                      char(track{i}.name)];

   nc.([name,'_longitude']).data=track{i}.lon;
   nc.([name,'_longitude']).units='degrees_east';
   nc.([name,'_longitude']).dimension={[name,'_length']};
   nc.([name,'_longitude']).long_name=['longitude, ',...
                                       char(track{i}.name)];

   nc.([name,'_distance']).data=track{i}.distance;
   nc.([name,'_distance']).units='nautical miles';
   nc.([name,'_distance']).dimension={[name,'_length']};
   nc.([name,'_distance']).long_name=['Distance along track ',...
                                      'from Bering Sea to ',...
                                      'Atlantic Ocean in 1NM spacing'];

   lats{i}=track{i}.lat;
   lons{i}=track{i}.lon;
   trackname{i}=char(track{i}.name);

  end

  ridgepack_write(nc,'InteRFACE_Shiptracks')

  shipshape=geoshape(lats,lons,'Name',trackname,'Geometry','line')

  shapewrite(shipshape,'InteRFACE_Shiptracks') 

  kmlwrite('InteRFACE_Shiptracks',shipshape,...
           'LineWidth',1,...
           'Color',[1 0.5 0],...
           'Name',trackname,...
           'Description','InteRFACE Ship Tracks at 1 NM spacing')

  ridgepack_fprint('png','InteRFACE_Shiptracks',1,2)

end

