
clear

laptop=true;

close all

if laptop

 gridlochr='/Volumes/MacBookProT3SeaIce/E3SM/highres/grid';
 datalochr='/Volumes/MacBookProT3SeaIce/E3SM/highres/monthly';

 gridlocv0='/Volumes/MacBookProT3SeaIce/E3SM/highresv0/grid';
 datalocv0='/Volumes/MacBookProT3SeaIce/E3SM/highresv0/monthly';

 gridloclr='/Volumes/MacBookProT3SeaIce/E3SM/lrhrequiv/grid';
 dataloclr='/Volumes/MacBookProT3SeaIce/E3SM/lrhrequiv/monthly';

 datalocobs='/Volumes/MacBookProT3SeaIce/';

 plotloc='/Volumes/MacBookProT3SeaIce/E3SM/highres/plots';

else

 gridlochr='/Volumes/CICESatArray/E3SM/highres/grid';
 datalochr='/Volumes/CICESatArray/E3SM/highres/monthly';

 gridlocv0='/Volumes/CICESatArray/E3SM/highresv0/grid';
 datalocv0='/Volumes/CICESatArray/E3SM/highresv0/monthly';

 gridloclr='/Volumes/CICESatArray/E3SM/lrhrequiv/grid';
 dataloclr='/Volumes/CICESatArray/E3SM/lrhrequiv/monthly';

 datalocobs='/Volumes/CICESatArray/data/SATELLITE/processed';

 plotloc='/Volumes/CICESatArray/E3SM/highres/plots';

end

alpha='abcdefghijklmnopqrstuvwxyz';

cont=[0 0.25 0.5 0.75 1.0 1.5:0.5:5 10];

k=0;
for rows=1:2
for cols=1:6

if cols<=3
 arctic=true;
else
 arctic=false;
end

% pull in grid information
if cols==1 | cols==4

 cd(gridlocv0)
 ncgrid=ridgepack_clone('hi.b1850c5_acmev0_highres.cice.h.0130-12.nc',...
                       {'TLAT','TLON','tarea'});

 cd(datalocv0)
 if rows==1
  nc=ridgepack_clone(...
       'b1850c5_acmev0_highres.HR_v0.state.mean.03.0031_0060.nc')
 elseif rows==2
  nc=ridgepack_clone(...
       'b1850c5_acmev0_highres.HR_v0.state.mean.09.0031_0060.nc')
 else
  error('wrong row')
 end

 nc=rmfield(nc,'hs');
 nc=rmfield(nc,'uvel');
 nc=rmfield(nc,'vvel');

 nc.volume=nc.hi;
 nc=rmfield(nc,'hi');

 nc.conc=nc.aice;
 nc=rmfield(nc,'aice');

 nc=ridgepack_struct(nc);

 nc.volume.data(nc.conc.data<0.001)=0;

 [cmap]=ridgepack_colormap(cont,0);

elseif cols==2 | cols==5

 cd(gridlochr)
 ncgrid=ridgepack_clone('gridFieldsRRS18to6v3.nc',{'latVertex','lonVertex',...
                        'verticesOnCell','indexToCellID','nEdgesOnCell'});
 nclat=ridgepack_clone('gridFieldsRRS18to6v3.nc',{'latCell'});
 nclon=ridgepack_clone('gridFieldsRRS18to6v3.nc',{'lonCell'});
 landm = shaperead('E3SM_HR_V1_C_grid_Coast.shp','UseGeoCoords',true);
 cd(datalochr)
 if rows==1
  nc=ridgepack_clone('mpascice.hist.am.HR_V1.mean.03.0026_0055.nc')
 elseif rows==2
  nc=ridgepack_clone('mpascice.hist.am.HR_V1.mean.09.0026_0055.nc')
 else
  error('wrong row')
 end

 nc.volume=nc.Monthly_avg_iceVolumeCell;
 nc=rmfield(nc,'Monthly_avg_iceVolumeCell');

 nc.conc=nc.Monthly_avg_iceAreaCell;
 nc=rmfield(nc,'Monthly_avg_iceAreaCell');

 nc=ridgepack_struct(nc);

 nc.volume.data(nc.conc.data<0.001)=0;

 [cmap]=ridgepack_colormap(cont,0);

elseif cols==3 | cols==6

  cd(gridloclr)
  ncgrid=ridgepack_clone('E3SM_LR_V1_grid',{'latVertex','lonVertex',...
                         'verticesOnCell','indexToCellID','nEdgesOnCell'});
  nclat=ridgepack_clone('E3SM_LR_V1_grid',{'latCell'});
  nclon=ridgepack_clone('E3SM_LR_V1_grid',{'lonCell'});
  landm = shaperead('E3SM_V1_C_grid_Coast.shp','UseGeoCoords',true);
  cd(dataloclr)
  if rows==1
   nc=ridgepack_clone('mpascice.hist.am.LR_V1.state.mean.03.0026_0055.nc')
  elseif rows==2
   nc=ridgepack_clone('mpascice.hist.am.LR_V1.state.mean.09.0026_0055.nc')
  else
   error('wrong row')
  end

  nc.volume=nc.timeMonthly_avg_iceVolumeCell;
  nc=rmfield(nc,'timeMonthly_avg_iceVolumeCell');

  nc.conc=nc.timeMonthly_avg_iceAreaCell;
  nc=rmfield(nc,'timeMonthly_avg_iceAreaCell');

  nc=ridgepack_struct(nc);

  nc.volume.data(nc.conc.data<0.001)=0;

  [cmap]=ridgepack_colormap(cont,0);

end

k=k+1;

ridgepack_multiplot(2,6,rows,cols,alpha(k));

if arctic

 if cols==1 & rows==1
  ridgepack_polarm('seaice','noland','grid','label')
 else
  ridgepack_polarm('seaice','noland')
 end

else

 if cols==4 & rows==1
  ridgepack_polarm('antarctic','noland','grid','label')
 else
  ridgepack_polarm('antarctic','noland')
 end

end

if cols==1 | cols==4

  for j=1:length(cont)-1

    c=[];
    d=[];
    e=[];

    if arctic
     if j<length(cont)
       idxn=find(nc.volume.data>cont(j) & ...
                 nc.volume.data<=cont(j+1) & ...
                 ncgrid.latitude.data>40);
     else
       idxn=find(nc.volume.data>cont(j+1) & ...
                 ncgrid.latitude.data>40);
     end
    else
     if j<length(cont)
       idxn=find(nc.volume.data>cont(j) & ...
                 nc.volume.data<=cont(j+1) & ...
                 ncgrid.latitude.data<-40);
     else
       idxn=find(nc.volume.data>cont(j+1) & ...
                 ncgrid.latitude.data<-40);
     end
    end

    [zindex,truecolor]=ridgepack_colorindex(nc.volume.data(idxn),cont,0);

    if length(idxn)>0

      maxidx=4;
      
      lat=zeros(length(idxn),maxidx+1);
      lon=zeros(length(idxn),maxidx+1);

      lats=nc.latt_bounds.data(1:maxidx,idxn)';
      lons=nc.lont_bounds.data(1:maxidx,idxn)';

      lat(:,1:maxidx)=lats;
      lon(:,1:maxidx)=lons;

      lat(:,maxidx+1)=lats(:,1);
      lon(:,maxidx+1)=lons(:,1);

      [cc,dd] = mfwdtran(gcm,lat(:,:),lon(:,:));

      ee = nc.conc.data(idxn)/100;

      patch(cc',dd',truecolor(1,:),'FaceVertexAlphaData',ee,...
                        'FaceAlpha','Flat','EdgeColor','none');

      drawnow

     end

     clear zindex truecolor cc dd ee lon lat idxn idyn

  end

  nc.mask=nc.volume;
  nc.mask.data(~isnan(nc.volume.data))=1;
  nc.mask.data(isnan(nc.volume.data))=0;

  ridgepack_maskm(nc.latitude.data,...
                  nc.longitude.data,...
                  nc.mask.data)

  drawnow

  if rows==1 & cols==1
   ylabel('March','FontSize',10)
  elseif rows==2 & cols==1
   ylabel('September','FontSize',10)
  end

  if rows==1 
   title('v0 HR','FontSize',10) 
  end

elseif cols==2 | cols==5

  for j=1:length(cont)-1

    if arctic
     if j<length(cont)
       idxn=find(nc.volume.data>cont(j) & ...
                 nc.volume.data<=cont(j+1) & ...
                 nclat.latCell.data>40*pi/180);
     else
       idxn=find(nc.volume.data>cont(j+1) & ...
                 nclat.latCell.data>40*pi/180);
     end
    else
     if j<length(cont)
       idxn=find(nc.volume.data>cont(j) & ...
                 nc.volume.data<=cont(j+1) & ...
                 nclat.latCell.data<-40*pi/180);
     else
       idxn=find(nc.volume.data>cont(j+1) & ...
                 nclat.latCell.data<-40*pi/180);
     end
    end

    [zindex,truecolor]=ridgepack_colorindex(nc.volume.data(idxn),cont,0);

    if length(idxn)>0

      lat=zeros(length(idxn),8);
      lon=zeros(length(idxn),8);

      for i=1:length(idxn)

       maxidx=ncgrid.nEdgesOnCell.data(idxn(i));

       lat(i,1:maxidx)=ncgrid.latitude.data(...
               ncgrid.verticesOnCell.data(...
                1:maxidx,idxn(i)))*180/pi;
       lon(i,1:maxidx)=ncgrid.longitude.data(...
               ncgrid.verticesOnCell.data(...
                1:maxidx,idxn(i)))*180/pi;
       
       lon(i,maxidx+1:8)=lon(i,1);
       lat(i,maxidx+1:8)=lat(i,1);

      end

      [cc,dd] = mfwdtran(gcm,lat(:,:),lon(:,:));

      ee = nc.conc.data(idxn);

      patch(cc',dd',truecolor(1,:),'FaceVertexAlphaData',ee,...
                         'FaceAlpha','Flat','EdgeColor','none')

      drawnow

     end

     clear zindex truecolor cc dd ee lon lat idxn idyn

  end

  [c,d] = mfwdtran(gcm,[landm.Lat],[landm.Lon]);
  h1=line(c,d,'Color',[0.5 0.5 0.5]);
  clear c d landm nc

  if rows==1
   title('v1 HR','FontSize',10) 
  end

elseif cols==3 | cols==6

  for j=1:length(cont)-1

    c=[];
    d=[];
    e=[];

    if arctic
     if j<length(cont)
       idxn=find(nc.volume.data>cont(j) & ...
                 nc.volume.data<=cont(j+1) & ...
                 nclat.latCell.data>40*pi/180);
     else
       idxn=find(nc.volume.data>cont(j+1) & ...
                 nclat.latCell.data>40*pi/180);
     end
    else
     if j<length(cont)
       idxn=find(nc.volume.data>cont(j) & ...
                 nc.volume.data<=cont(j+1) & ...
                 nclat.latCell.data<-40*pi/180);
     else
       idxn=find(nc.volume.data>cont(j+1) & ...
                 nclat.latCell.data<-40*pi/180);
     end
    end

    [zindex,truecolor]=ridgepack_colorindex(nc.volume.data(idxn),cont,0);

    if length(idxn)>0

      lat=zeros(length(idxn),8);
      lon=zeros(length(idxn),8);

      for i=1:length(idxn)

       maxidx=ncgrid.nEdgesOnCell.data(idxn(i));

       lat(i,1:maxidx)=ncgrid.latitude.data(...
               ncgrid.verticesOnCell.data(...
                1:maxidx,idxn(i)))*180/pi;
       lon(i,1:maxidx)=ncgrid.longitude.data(...
               ncgrid.verticesOnCell.data(...
                1:maxidx,idxn(i)))*180/pi;

       lon(i,maxidx+1:8)=lon(i,1);
       lat(i,maxidx+1:8)=lat(i,1);

      end

      [cc,dd] = mfwdtran(gcm,lat(:,:),lon(:,:));

      ee = nc.conc.data(idxn);

      patch(cc',dd',truecolor(1,:),'FaceVertexAlphaData',ee,...
                       'FaceAlpha','Flat','EdgeColor','none')

      drawnow

     end

     clear zindex truecolor cc dd ee lon lat idxn idyn

  end

  [c,d] = mfwdtran(gcm,[landm.Lat],[landm.Lon]);
  h1=line(c,d,'Color',[0.5 0.5 0.5]);
  clear c d landm nc

  if rows==1
   title('v1 LR','FontSize',10) 
  end

  if rows==1 & cols==3
   ridgepack_colorbar(cont,'m')
   ridgepack_cbshare(gca)
  end

end


% plot ice edge in orange
cd(datalocobs)

if rows==1
 if arctic
  ncconc=ridgepack_clone('G02202_v3_merged_conc_north_1979_1999_Mar_mean');
 else
  ncconc=ridgepack_clone('G02202_v3_merged_conc_south_1979_1999_Mar_mean');
 end
elseif rows==2
 if arctic
  ncconc=ridgepack_clone('G02202_v3_merged_conc_north_1979_1999_Sep_mean');
 else
  ncconc=ridgepack_clone('G02202_v3_merged_conc_south_1979_1999_Sep_mean');
 end
else
 error('rows')
end

conc=ncconc.conc.data+ncconc.conc_std.data;
conc(conc<15)=0;
conc(conc>14)=1;
conc(ncconc.latitude.data>=84)=1;
h1=ridgepack_maskm(ncconc.latitude.data,ncconc.longitude.data,...
                      conc,[0.9100 0.4100 0.1700],0.2);

conc=ncconc.conc.data-ncconc.conc_std.data;
conc(conc<15)=0;
conc(conc>14)=1;
conc(ncconc.latitude.data>=84)=1;
h1=ridgepack_maskm(ncconc.latitude.data,ncconc.longitude.data,...
                      conc,[0.9100 0.4100 0.1700],0.2);

conc=ncconc.conc.data;
conc(conc<15)=0;
conc(conc>14)=1;
conc(ncconc.latitude.data>=84)=1;
h1=ridgepack_maskm(ncconc.latitude.data,ncconc.longitude.data,...
                      conc,[0.9100 0.4100 0.1700]);

drawnow

end
end

ridgepack_multialign(gcf)

cd(plotloc)

ridgepack_fprint('png','E3SM_HR_V1_Sea_Ice_Thickness',1,2);
ridgepack_fprint('epsc','E3SM_HR_V1_Sea_Ice_Thickness',1,2);


