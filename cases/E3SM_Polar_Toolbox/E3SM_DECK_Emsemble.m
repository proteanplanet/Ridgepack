
clear
clf


% pull in grid information
cd('/Users/aroberts/work/E3SM')
ncgrid=ridgepack_clone('oEC60to30v3_60layer.restartFrom_anvil0926.171101',...
                       {'latVertex','lonVertex','verticesOnCell',...
                        'indexToCellID','nEdgesOnCell'});

nclat=ridgepack_clone('oEC60to30v3_60layer.restartFrom_anvil0926.171101',...
                       {'latCell','areaCell'});

notation={'(a)','(c)','(b)','(d)'};

% hemisphere
for hem=1:2

 clf

 timebrackets=[1 2];

 place=0;

 for timebracket=timebrackets


 if timebracket==1
  startyear=1;
  endyear=20;  
 elseif timebracket==2
  startyear=21;
  endyear=35;
 end

 spanyear=abs(endyear-startyear)+1;

 if hem==1
  filename='E3SM_V1_North'
 elseif hem==2
  filename='E3SM_V1_South'
 end

 if startyear==1 & endyear==20
  header='1980-1999';
  yearstart=1980;
  yearend=1999;
 elseif startyear==21 & endyear==35
  header='2000-2014';
  yearstart=2000;
  yearend=2014;
 end

 if hem==1
  ncfile='G02202_v3_merged_conc_north_1979_2017_time_bounds.nc';
 elseif hem==2
  ncfile='G02202_v3_merged_conc_south_1979_2017_time_bounds.nc';
 end
 fieldoa='conc';

 ensembles=[1 2 3 4 5];

 for month=1:2

  place=place+1;

  if month==1
   cd('/Users/aroberts/data/SATELLITE/processed')
   [ncconc]=ridgepack_timesubset(ncfile,'conc',[3],yearstart,yearend)
   cd('/Users/aroberts/work/andrewIceData/marchIce')
  elseif month==2
   cd('/Users/aroberts/data/SATELLITE/processed')
   [ncconc]=ridgepack_timesubset(ncfile,'conc',[9],yearstart,yearend)
   cd('/Users/aroberts/work/andrewIceData/septIce')
  end

  ncconc.(fieldoa).data(ncconc.(fieldoa).data<15)=0;
  ncconc.(fieldoa).data(ncconc.(fieldoa).data>14)=1;
  if startyear<1988
   ncconc.(fieldoa).data(ncconc.latitude.data>=84)=1;
  else
   ncconc.(fieldoa).data(ncconc.latitude.data>87)=1;
  end

  for ensemblemember=ensembles;

   for i=startyear:endyear
    ncx=ridgepack_reduce(...
      ridgepack_clone(['historical',num2str(ensemblemember)],...
       {'timeMonthly_avg_iceVolumeCell',...
          'timeMonthly_avg_iceAreaCell'},i),{'time'});
    if i==startyear
     nc=ncx;
    else
     nc.timeMonthly_avg_iceVolumeCell.data=...
      nc.timeMonthly_avg_iceVolumeCell.data+...
       ncx.timeMonthly_avg_iceVolumeCell.data
     nc.timeMonthly_avg_iceAreaCell.data=...
      nc.timeMonthly_avg_iceAreaCell.data+...
       ncx.timeMonthly_avg_iceAreaCell.data
    end
   end

   nc.timeMonthly_avg_iceVolumeCell.data=...
    nc.timeMonthly_avg_iceVolumeCell.data./spanyear;

   nc.timeMonthly_avg_iceAreaCell.data=...
    nc.timeMonthly_avg_iceAreaCell.data./spanyear;

   % calculate area weighted thickness
   if hem==1
    idxvol=find(nc.timeMonthly_avg_iceVolumeCell.data>0 & ...
                nc.timeMonthly_avg_iceAreaCell.data>0.15 & ...
                nclat.latCell.data>0);
   elseif hem==2
    idxvol=find(nc.timeMonthly_avg_iceVolumeCell.data>0 & ...
                nc.timeMonthly_avg_iceAreaCell.data>0.15 & ...
                nclat.latCell.data<0);
   end
   totalweight=sum(nclat.areaCell.data(idxvol));
   vol=sum(nc.timeMonthly_avg_iceVolumeCell.data(idxvol).* ...
           nclat.areaCell.data(idxvol))./totalweight;

   thickness(ensemblemember)=vol;

   if ensemblemember==ensembles(1)
    ncv=nc;

    ncv=rmfield(nc,'attributes');
    ncv.attributes.title='E3SM';

    ncv.volume=nc.timeMonthly_avg_iceVolumeCell;

    ncv.conc=nc.timeMonthly_avg_iceAreaCell;
   else
    ncv.volume.data=ncv.volume.data+nc.timeMonthly_avg_iceVolumeCell.data;

    ncv.conc.data=ncv.conc.data+nc.timeMonthly_avg_iceAreaCell.data;
   end

   clear ncx nc;

  end
   
  ncv=ridgepack_struct(ncv);

  ncv.volume.data=ncv.volume.data/length(ensembles);

  ncv.conc.data=ncv.conc.data/length(ensembles);

  ncv.volume.data(ncv.conc.data<0.15)=0;

  % calculate area weighted thickness
  if hem==1
   idxvol=find(ncv.volume.data>0 & nclat.latCell.data>0);
  elseif hem==2
   idxvol=find(ncv.volume.data>0 & nclat.latCell.data<0);
  end
  totalweight=sum(nclat.areaCell.data(idxvol));
  vol=sum(ncv.volume.data(idxvol).*nclat.areaCell.data(idxvol))./totalweight;

  thicknessmean=vol;

  cont=[0:0.5:5];

  [cmap]=ridgepack_colormap(cont,0);

  ridgepack_multiplot(2,length(timebrackets),month,timebracket,...
                        char(notation{place}))

  if timebracket==1 & month==1
   if hem==1
    ridgepack_polarm('seaice','grid','label')
   elseif hem==2
    ridgepack_polarm('antarctic','grid','label')
   end
  else
   if hem==1
    ridgepack_polarm('seaice','grid')
   elseif hem==2
    ridgepack_polarm('antarctic','grid')
   end
  end

  mstruct=gcm;

  for j=1:length(cont)-1

   if hem==1
     if j<length(cont)
      idxn=find(ncv.volume.data>cont(j) & ncv.volume.data<=cont(j+1) & ...
                nclat.latCell.data>0);
     else
      idxn=find(ncv.volume.data>cont(j+1) & nclat.latCell.data>0);
     end
   elseif hem==2
     if j<length(cont)
      idxn=find(ncv.volume.data>cont(j) & ncv.volume.data<=cont(j+1) & ...
                nclat.latCell.data<0);
     else
      idxn=find(ncv.volume.data>cont(j+1) & nclat.latCell.data<0);
     end
   else
     error('hemisphere error')
   end

   [zindex,truecolor]=ridgepack_colorindex(ncv.volume.data(idxn),cont,0);
   if length(idxn)>0

     lat=zeros(length(idxn),8);
     lon=zeros(length(idxn),8);

     px=[];
     py=[];

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

      [c,d] = mfwdtran(mstruct,lat(i,:),lon(i,:),gca,'surface');
      
      px=[px; NaN; c];
      py=[py; NaN; d];

     end

     plot(px,py,'Color',truecolor(1,:),'LineWidth',0.5)

     drawnow

    end

   end

   h1=ridgepack_maskm(ncconc.latitude.data,ncconc.longitude.data,...
                      ncconc.(fieldoa).data,'m',0.4);
 
   if timebracket==1 & month==1
    [hcb]=ridgepack_colorbar(cont,'m')
    ridgepack_multicb
   end

   if month==1
    title([header],'FontSize',11)
   end

   if timebracket==1 & month==1
    ylabel('March','FontSize',11)
   elseif timebracket==1 & month==2
    ylabel('September','FontSize',11)
   end

   h1=text(max(xlim),min(ylim),[num2str(min(thickness),'%4.2f'),' m '],...
        'FontSize',10,... 
        'Color','b',...
        'Margin',1,...
        'BackgroundColor','w',...
        'HorizontalAlignment','right',...
        'VerticalAlignment','bottom')  

   ext1=get(h1,'extent');

   h2=text(max(xlim),min(ylim)+ext1(4),...
        [num2str(thicknessmean,'%4.2f'),' m '],...
        'FontSize',10,... 
        'Color','k',...
        'Margin',1,...
        'BackgroundColor','w',...
        'HorizontalAlignment','right',...
        'VerticalAlignment','bottom')  

   ext2=get(h2,'extent');

   text(max(xlim),min(ylim)+ext1(4)+ext2(4),...
        [num2str(max(thickness),'%4.2f'),' m '],...
        'FontSize',10,... 
        'Color','r',...
        'Margin',1,...
        'BackgroundColor','w',...
        'HorizontalAlignment','right',...
        'VerticalAlignment','bottom')  

   drawnow

   disp('-----------------------------------------------')
   disp('Mean NH Thickness ')
   thickness
   thicknessmean
   disp('-----------------------------------------------')
 

  end

 end

 ridgepack_multialign(gcf);

 cd('/Users/aroberts/work/andrewIceData')

 ridgepack_fprint('png',filename,1,2)
 ridgepack_fprint('epsc',filename,1,2)

end



