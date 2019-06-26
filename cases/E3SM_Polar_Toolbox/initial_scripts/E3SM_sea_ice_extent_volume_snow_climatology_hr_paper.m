
close all
clear

laptop=true;
%laptop=false;

colss=colormap(lines(4));
cols=colss;
cols(1,:)=colss(2,:);
cols(2,:)=colss(3,:);
cols(3,:)=colss(1,:);
cols(4,:)=colss(4,:);

list={'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'};

alpha='abcdefghijklmnopqrstuvwxyz';

k=0

for col=1:2

if laptop
 cd('/Volumes/MacBookProT3SeaIce/E3SM/highres/monthly')
else
 cd('/Users/afroberts/SIhMSatArray/E3SM/highres/monthly')
end
if col==1
 nchrn=ridgepack_clone('E3SM-HR-V1_mpascice.hist.am.timeSeriesStatsMonthly.north_volume_area_extent')
else
 nchrn=ridgepack_clone('E3SM-HR-V1_mpascice.hist.am.timeSeriesStatsMonthly.south_volume_area_extent')
end

if laptop
 cd('/Volumes/MacBookProT3SeaIce/E3SM/lrhrequiv/monthly')
else
 cd('/Users/afroberts/SIhMSatArray/E3SM/lrhrequiv/monthly')
end
if col==1
 nclrn=ridgepack_clone('E3SM-LR-V1_mpascice.hist.am.timeSeriesStatsMonthly.north_volume_area_extent')
else
 nclrn=ridgepack_clone('E3SM-LR-V1_mpascice.hist.am.timeSeriesStatsMonthly.south_volume_area_extent')
end

if laptop
 cd('/Volumes/MacBookProT3SeaIce/E3SM/highresv0/monthly')
else
 cd('/Users/afroberts/SIhMSatArray/E3SM/v0highres/monthly')
end
if col==1
 nchrnv0=ridgepack_clone('b1850c5_acmev0_highres.cice.h.Northern_volume_area_extent.nc')
else
 nchrnv0=ridgepack_clone('b1850c5_acmev0_highres.cice.h.Southern_volume_area_extent.nc')
end

if laptop
 cd('/Volumes/MacBookProT3SeaIce')
else
 cd('/Users/afroberts/SIhMSatArray/data/SATELLITE/processed/')
end
if col==1
 ncobs=ridgepack_clone('G02202_v3_merged_conc_north_1979_1999_monthly_mean')
else
 ncobs=ridgepack_clone('G02202_v3_merged_conc_south_1979_1999_monthly_mean')
end



dvhr=datevec(nchrn.time.data);
yearshr=dvhr(:,1);
idxn=find(yearshr>=26 & yearshr<=55);
timehr=nchrn.time.data(idxn);
extenthr=nchrn.extent.data(idxn)/10^6;
volumehr=nchrn.volume.data(idxn)/10^3;
snowhr=nchrn.snowvolume.data(idxn)/10^2;
dvhr=datevec(timehr);
monthshr=dvhr(:,2);

dvv0=datevec(nchrnv0.time.data);
yearsv0=dvv0(:,1);
idxn=find(yearsv0>=31 & yearsv0<=60);
timev0=nchrnv0.time.data(idxn);
extentv0=nchrnv0.extent.data(idxn)/10^6;
volumev0=nchrnv0.volume.data(idxn)/10^3;
snowv0=nchrnv0.snowvolume.data(idxn)/10^2;
dvv0=datevec(timev0);
monthsv0=dvv0(:,2);

dvlr=datevec(nclrn.time.data);
yearslr=dvlr(:,1);
idxn=find(yearslr>=26 & yearslr<=55);
timelr=nclrn.time.data(idxn);
extentlr=nclrn.extent.data(idxn)/10^6;
volumelr=nclrn.volume.data(idxn)/10^3;
snowlr=nclrn.snowvolume.data(idxn)/10^2;
dvlr=datevec(timelr);
monthslr=dvlr(:,2);

dvob=datevec(ncobs.time.data);
yearsob=dvob(:,1);
idxn=find(yearsob>=1979 & yearsob<=1999);
timeob=ncobs.time.data(idxn);
extentob=ncobs.extent.data(idxn)/10^6;
dvob=datevec(timeob);
monthsob=dvob(:,2);

k=k+1;

ridgepack_multiplot(3,2,1,col,alpha(k))
for i=1:12

 idxn=find(monthsv0==i);
 medev0(i)=median(extentv0(idxn));
 maxev0(i)=max(extentv0(idxn));
 minev0(i)=min(extentv0(idxn));
 h0=plot([i-0.15 i-0.15],[minev0(i) maxev0(i)],'Color',cols(3,:));
 hold on
 plot(i-0.15,medev0(i),'.','Color',cols(3,:),'MarkerSize',7.5)
 plot([i-0.15-0.05 i-0.15+0.05],[maxev0(i) maxev0(i)],'Color',cols(3,:))
 plot([i-0.15-0.05 i-0.15+0.05],[minev0(i) minev0(i)],'Color',cols(3,:))

 idxn=find(monthslr==i);
 medelr(i)=median(extentlr(idxn));
 maxelr(i)=max(extentlr(idxn));
 minelr(i)=min(extentlr(idxn));
 h3=plot([i+0.15 i+0.15],[minelr(i) maxelr(i)],'Color',cols(2,:));
 hold on
 plot(i+0.15,medelr(i),'.','Color',cols(2,:),'MarkerSize',7.5)
 plot([i+0.15-0.05 i+0.15+0.05],[maxelr(i) maxelr(i)],'Color',cols(2,:))
 plot([i+0.15-0.05 i+0.15+0.05],[minelr(i) minelr(i)],'Color',cols(2,:))

 idxn=find(monthsob==i);
 medeob(i)=median(extentob(idxn));
 maxeob(i)=max(extentob(idxn));
 mineob(i)=min(extentob(idxn));
 h2=plot([i+0.05 i+0.05],[mineob(i) maxeob(i)],'Color',cols(4,:));
 hold on
 plot(i+0.05,medeob(i),'.','Color',cols(4,:),'MarkerSize',10.5)
 plot([i+0.05-0.05 i+0.05+0.05],[maxeob(i) maxeob(i)],'Color',cols(4,:))
 plot([i+0.05-0.05 i+0.05+0.05],[mineob(i) mineob(i)],'Color',cols(4,:))

 idxn=find(monthshr==i);
 medehr(i)=median(extenthr(idxn));
 maxehr(i)=max(extenthr(idxn));
 minehr(i)=min(extenthr(idxn));
 h1=plot([i-0.05 i-0.05],[minehr(i) maxehr(i)],'Color',cols(1,:));
 hold on
 plot(i-0.05,medehr(i),'.','Color',cols(1,:),...
         'MarkerSize',10.5,'MarkerFaceColor',cols(1,:));
 plot([i-0.05-0.05 i-0.05+0.05],[maxehr(i) maxehr(i)],'Color',cols(1,:))
 plot([i-0.05-0.05 i-0.05+0.05],[minehr(i) minehr(i)],'Color',cols(1,:))

end

xx=[1:0.1:12];
yy=spline([1:12]',medeob',xx);
hh=plot(xx,yy,':')
set(hh,'Color',cols(4,:))

yy=spline([1:12]',medehr',xx);
hh=plot(xx,yy,':')
set(hh,'Color',cols(1,:))

%plot([1:12],medehr,'Color',cols(1,:),'LineStyle',':')
%plot([1:12],medelr,'Color',cols(2,:),'LineStyle',':')
%plot([1:12],medeob,'Color',cols(3,:),'LineStyle',':')

if col==2
 legend([h0 h1 h3 h2],{'v0 HR','v1 HR','v1 LR','NOAA CDR'},'Location','SouthEast')
 legend('boxoff')
else
 legend('off')
end

ylim([0 25])
xlim([0.5 12.5])
set(gca,'Xtick',[1:1:12],'XTickLabel',[])
if col==2
 set(gca,'YTickLabel',[])
 title('Southern Ocean Sea Ice')
 ylabel('')
 xlabel('')
else
 ylabel('$\mathrm{Extent\,(\times\,10^{6}\,km^{2})}$')
 title('Arctic System Sea Ice')
 xlabel('')
end

k=k+1;

ridgepack_multiplot(3,2,2,col,alpha(k))
for i=1:12

 idxn=find(monthsv0==i);
 medvv0(i)=median(volumev0(idxn));
 maxvv0(i)=max(volumev0(idxn));
 minvv0(i)=min(volumev0(idxn));
 plot([i-0.1 i-0.1],[minvv0(i) maxvv0(i)],'Color',cols(3,:))
 hold on
 plot(i-0.1,medvv0(i),'.','Color',cols(3,:),'MarkerSize',7.5)
 plot([i-0.1-0.05 i-0.1+0.05],[maxvv0(i) maxvv0(i)],'Color',cols(3,:))
 plot([i-0.1-0.05 i-0.1+0.05],[minvv0(i) minvv0(i)],'Color',cols(3,:))

 idxn=find(monthslr==i);
 medvlr(i)=median(volumelr(idxn));
 maxvlr(i)=max(volumelr(idxn));
 minvlr(i)=min(volumelr(idxn));
 plot([i+0.1 i+0.1],[minvlr(i) maxvlr(i)],'Color',cols(2,:))
 hold on
 plot(i+0.1,medvlr(i),'.','Color',cols(2,:),'MarkerSize',7.5)
 plot([i+0.1-0.05 i+0.1+0.05],[maxvlr(i) maxvlr(i)],'Color',cols(2,:))
 plot([i+0.1-0.05 i+0.1+0.05],[minvlr(i) minvlr(i)],'Color',cols(2,:))

 idxn=find(monthshr==i);
 medvhr(i)=median(volumehr(idxn));
 maxvhr(i)=max(volumehr(idxn));
 minvhr(i)=min(volumehr(idxn));
 plot([i i],[minvhr(i) maxvhr(i)],'Color',cols(1,:))
 hold on
 plot(i,medvhr(i),'.','Color',cols(1,:),'MarkerSize',10.5)
 plot([i-0.05 i+0.05],[maxvhr(i) maxvhr(i)],'Color',cols(1,:))
 plot([i-0.05 i+0.05],[minvhr(i) minvhr(i)],'Color',cols(1,:))

end

yy=spline([1:12]',medvhr',xx);
hh=plot(xx,yy,':')
set(hh,'Color',cols(1,:))

ylim([0 46])
xlim([0.5 12.5])
set(gca,'Xtick',[1:1:12],'XTickLabel',[])
if col==2
 set(gca,'YTickLabel',[])
 ylabel('')
 xlabel('')
else
 ylabel('$\mathrm{Volume\,(\times\,10^{3}\,km^{3})}$')
 xlabel('')
end

k=k+1;

ridgepack_multiplot(3,2,3,col,alpha(k))
for i=1:12

 idxn=find(monthsv0==i);
 medsv0(i)=median(snowv0(idxn));
 maxsv0(i)=max(snowv0(idxn));
 minsv0(i)=min(snowv0(idxn));
 plot([i-0.1 i-0.1],[minsv0(i) maxsv0(i)],'Color',cols(3,:))
 hold on
 plot(i-0.1,medsv0(i),'.','Color',cols(3,:),'MarkerSize',7.5)
 plot([i-0.1-0.05 i-0.1+0.05],[maxsv0(i) maxsv0(i)],'Color',cols(3,:))
 plot([i-0.1-0.05 i-0.1+0.05],[minsv0(i) minsv0(i)],'Color',cols(3,:))

 idxn=find(monthslr==i);
 medslr(i)=median(snowlr(idxn));
 maxslr(i)=max(snowlr(idxn));
 minslr(i)=min(snowlr(idxn));
 plot([i+0.1 i+0.1],[minslr(i) maxslr(i)],'Color',cols(2,:))
 hold on
 plot(i+0.1,medslr(i),'.','Color',cols(2,:),'MarkerSize',7.5)
 plot([i+0.1-0.05 i+0.1+0.05],[maxslr(i) maxslr(i)],'Color',cols(2,:))
 plot([i+0.1-0.05 i+0.1+0.05],[minslr(i) minslr(i)],'Color',cols(2,:))

 idxn=find(monthshr==i);
 medshr(i)=median(snowhr(idxn));
 maxshr(i)=max(snowhr(idxn));
 minshr(i)=min(snowhr(idxn));
 plot([i i],[minshr(i) maxshr(i)],'Color',cols(1,:))
 hold on
 plot(i,medshr(i),'.','Color',cols(1,:),'MarkerSize',10.5)
 plot([i-0.05 i+0.05],[maxshr(i) maxshr(i)],'Color',cols(1,:))
 plot([i-0.05 i+0.05],[minshr(i) minshr(i)],'Color',cols(1,:))

end

yy=spline([1:12]',medshr',xx);
hh=plot(xx,yy,':')
set(hh,'Color',cols(1,:))

ylim([-2 69])
xlim([0.5 12.5])
set(gca,'Xtick',[1:1:12])%,'XTickLabel',list)
if col==2
 set(gca,'YTickLabel',[])
 ylabel('')
else
 ylabel('$\mathrm{Snow Volume\,(\times\,10^{2}\,km^{3})}$')
end
xlabel('Month')

end % col

ridgepack_multialign(gcf)

if laptop
 cd('/Volumes/MacBookProT3SeaIce/E3SM/highres/plots')
else
 cd('/Users/afroberts/SIhMSatArray/E3SM/highres/plots') 
end

ridgepack_fprint('png','sea_ice_climate',1,2)
ridgepack_fprint('epsc','sea_ice_climate',1,2)






