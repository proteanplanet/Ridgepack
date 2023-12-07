

clear
close all

filetabe='Icedge';

legnames={'PSLV','Icedge'};

clear lr na

cd(['/Users/afroberts/data/MODEL/E3SM/PSLV']);
nc1=ridgepack_clone('PSLV.ensemble.lemnisc.0051-0100.PI.nc','totalKineticEnergy');
divisions=ceil(max(nc1.totalKineticEnergy.data(1,:))./200);
xbins=[0:divisions:divisions*200]./10^2;
h1=histogram(nc1.totalKineticEnergy.data(1,:)./10^2,xbins,'Normalization','probability')
xbins=h1.BinEdges;
lr{1}=h1.Values;
eta(1)=median(nc1.totalKineticEnergy.data(1,:))/10^2;
sum(h1.Values)
h2=histogram(nc1.totalKineticEnergy.data(2,:)./10^2,xbins,'Normalization','probability')
lr{2}=h2.Values;
eta(2)=median(nc1.totalKineticEnergy.data(2,:))./10^2;
sum(h2.Values)
h3=histogram(nc1.totalKineticEnergy.data(3,:)./10^2,xbins,'Normalization','probability')
lr{3}=h3.Values;
eta(3)=median(nc1.totalKineticEnergy.data(3,:))./10^2;
sum(h3.Values)
mrp=max(nc1.totalKineticEnergy.data(3,:))./10^2

clf

cd(['/Users/afroberts/data/MODEL/E3SM/Icedge']);
nc1=ridgepack_clone('Icedge.ensemble.lemnisc.0051-0100.PI.nc','totalKineticEnergy');
h1=histogram(nc1.totalKineticEnergy.data(1,:)./10^2,xbins,'Normalization','probability')
xbins=h1.BinEdges;
lr{4}=h1.Values;
eta(4)=median(nc1.totalKineticEnergy.data(1,:))./10^2;
sum(h1.Values)
h2=histogram(nc1.totalKineticEnergy.data(2,:)./10^2,xbins,'Normalization','probability')
lr{5}=h2.Values;
eta(5)=median(nc1.totalKineticEnergy.data(2,:))./10^2;
sum(h2.Values)
h3=histogram(nc1.totalKineticEnergy.data(3,:)./10^2,xbins,'Normalization','probability')
lr{6}=h3.Values;
eta(6)=median(nc1.totalKineticEnergy.data(3,:))./10^2;
sum(h3.Values)
mrl=max(nc1.totalKineticEnergy.data(3,:))./10^2


clf

ahs='abcd';

cols=lines(3);
for colu=1:2

 ridgepack_multiplot(1,2,1,colu,ahs(colu));
 if colu==1
  traces=[2 3 5 6];
  yticks=[0:0.01:0.08];
 else
  traces=[1 4];
  yticks=[0:0.01:0.03];
 end

 for k=traces
  if k==1 | k==4
   col=cols(1,:);
  elseif k==2 | k==5
   col=cols(2,:);
  elseif k==3 | k==6
   col=cols(3,:);
  end
  clear stepx stepy
  for j=1:length(xbins)
   if j==1;
    stepx(1)=xbins(j);
    stepy(1)=0;
    stepx(end+1)=xbins(j);
    stepy(end+1)=lr{k}(j);
   elseif j==length(xbins);
    stepx(end+1)=xbins(j);
    stepy(end+1)=lr{k}(j-1);
    stepx(end+1)=xbins(j);
    stepy(end+1)=0;
   else
    stepx(end+1)=xbins(j);
    stepy(end+1)=lr{k}(j-1);
    stepx(end+1)=xbins(j);
    stepy(end+1)=lr{k}(j);
   end
  end
  if k<4
   h(k)=plot(stepx,stepy,':','Color',col)
  else
   h(k)=plot(stepx,stepy,'-','Color',col)
  end
  hold on
 end

 if colu==1
  xlim([0 10])
  ylim([0 0.065])
  set(gca,'Ytick',yticks)
  legend([h(6) h(3) h(5) h(2)],...
         {['\eta=',num2str(eta(6),'%2.2f'),' SH Icedge'],...
          ['\eta=',num2str(eta(3),'%2.2f'),' SH PSLV'],...
          ['\eta=',num2str(eta(5),'%2.2f'),' NH Icedge'],...
          ['\eta=',num2str(eta(2),'%2.2f'),' NH PSLV']})
  legend('boxoff')
  ylabel('Probability')
  xlabel('Total Sea Ice Kinetic Energy \times{10^2} TJ') 
  title('Hemispheric','FontWeight','normal')
 else
  xlim([0 10])
  ylim([0 0.025])
  set(gca,'Ytick',yticks)
  legend([h(4) h(1)],...
         {['\eta=',num2str(eta(4),'%2.2f'),', max=',num2str(mrl,'%2.2f'),' Icedge'],...
          ['\eta=',num2str(eta(1),'%2.2f'),', max=',num2str(mrp,'%2.2f'),' PSLV']})
  legend('boxoff')
  ylabel(' ')
  xlabel('Total Sea Ice Kinetic Energy \times{10^2} TJ') 
  title('Global','FontWeight','normal')
 end

end

ridgepack_multialign(gcf)

cd('/Users/afroberts/work')

ridgepack_fprint('png','PSLV_Icedge_KE_histogram',1,2)
