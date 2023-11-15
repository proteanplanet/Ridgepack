
% This script provides an example of rendering E3SM sea ice model data
% on a map in Ridgepack
% 
% Andrew Roberts, LANL, 2019 (afroberts@lanl.gov)

clear
close all

% location of model E3SM data   
dire3smdata=['/Volumes/MacBookProT3SeaIce/E3SM/DECK/monthly/'];

% location if E3SM grid data
gridloc='/Users/afroberts/data/E3SM/DECK/grid';

% print location
printloc=dire3smdata;

% set start months
springstartmonth=[3,3,3,3,3]
springendmonth=[5,5,5,5,5]

% set end months
fallstartmonth=[10,10,10,10,10];
fallendmonth=[12,12,12,12,12];

% set years
year=[2004:2008];

% alphabetical 
alpha='abcdefghijklmnopqrstuvwzyz';
alphak=0;


% grab data for season (1=spring, 2=fall)
for si=1:2

 % grab observed data
 if si==1
  startmonth=springstartmonth;
  endmonth=springendmonth;
 else
  startmonth=fallstartmonth;
  endmonth=fallendmonth;
 end

 % import base model grid data for E3SM
 cd(gridloc)
 ncvert=ridgepack_clone('E3SM_LR_V1_grid',{'latVertex',...
                                           'lonVertex',...
                                           'verticesOnCell',...
                                           'indexToCellID',...
                                           'nEdgesOnCell',...
                                           'edgesOnCell',...
                                           'cellsOnEdge'});

 nclat=ridgepack_clone('E3SM_LR_V1_grid',{'latCell'});

 mask=nclat.nCells.data(rad2deg(nclat.latCell.data)>40);

 % build seasonal ice thickness mean from E3SM over central arctic
 % for five ensemble members hi takes 1, 2, 3, 4 or 5
 maxensemble=3;
 for hi=1:maxensemble

  % construct mean from model data 
  count=0;
  for yi=1:length(year)
  for mi=startmonth(yi):endmonth(yi)
   count=count+1;
   filem=[dire3smdata,'h',num2str(hi),'/archive/ice/hist/',...
          'mpascice.hist.am.timeSeriesStatsMonthly.',...
          num2str(year(yi),'%4.4i'),'-',...
          num2str(mi,'%2.2i'),'-01.nc'];
   nc=ridgepack_clone(filem,{'timeMonthly_avg_iceVolumeCell',...
                             'timeMonthly_avg_iceAreaCell'});
   nc=ridgepack_reduce(nc,{'time'});
   if count==1
    nc.attributes=[];
    nc.attributes.title=['Ensemble ',num2str(hi),' E3SM data'];
    ncm=nc;
   else
    ncm.timeMonthly_avg_iceVolumeCell.data=...
             nc.timeMonthly_avg_iceVolumeCell.data + ...
             ncm.timeMonthly_avg_iceVolumeCell.data;
    ncm.timeMonthly_avg_iceAreaCell.data=...
             nc.timeMonthly_avg_iceAreaCell.data + ...
             ncm.timeMonthly_avg_iceAreaCell.data;
   end

  end
  end


  ncm.timeMonthly_avg_iceVolumeCell.data=...
              ncm.timeMonthly_avg_iceVolumeCell.data/count;
  ncm.timeMonthly_avg_iceAreaCell.data=...
              ncm.timeMonthly_avg_iceAreaCell.data/count;

  alphak=alphak+1;
  ridgepack_multiplot(2,maxensemble,si,hi,char(alpha(alphak)));
  ridgepack_polarm('centralarctic2','grid','label','noland')
  cont=[0:0.25:0.75 1.0:0.5:5.5];
  ref=0;
  if alphak>1
   ridgepack_pcole3sm(ncm,'timeMonthly_avg_iceVolumeCell',ncvert,mask,...
                      cont,'linear',ref,'none')
  else
   ridgepack_pcole3sm(ncm,'timeMonthly_avg_iceVolumeCell',ncvert,mask,...
                      cont,'linear',ref,'horizontal')
   ridgepack_multicb(gca)
  end
  ridgepack_pcoaste3sm(ncvert)
  if si==1
   title(['Ensemble ',num2str(hi)],'FontSize',10)
  else 
   title('')
  end
  if si==1 & hi==1
   ylabel('Spring','FontSize',10)
  elseif si==2 & hi==1
   ylabel('Fall','FontSize',10)
  end

 end

end

ridgepack_multialign(gcf,['E3SM Sea Ice Volume ',num2str(year(1)),'-',num2str(year(end))],...
                     11,[0 0 0],3)

cd(printloc)

ridgepack_fprint('png','test',1,1)

