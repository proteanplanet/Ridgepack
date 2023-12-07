
close all
clear

% plot location
plotloc='/Users/afroberts/work';
 
gridloc1=['/Users/afroberts/data/MODEL/E3SM/v3/IcoswISC30E3r2'];
gridfile1='mpaso.IcoswISC30E3r2.20230901.nc';
leg1='Unsmoothed Icos 30';

gridloc2=['/Users/afroberts/data/MODEL/E3SM/v3/IcoswISC30E3r6'];
gridfile2='initial_state.nc';
leg2='Smoothed Icos 30 at 100km';
  
% obtain grid information
cd(gridloc1)
nccell1=ridgepack_clone(gridfile1,{'latCell','lonCell','bottomDepth'});
  
cd(gridloc2)
nccell2=ridgepack_clone(gridfile2,{'latCell','lonCell','bottomDepth'});

latmasknorth=find(nccell1.latitude.data>60*pi/180);
latmasksouth=find(nccell1.latitude.data<-60*pi/180);

% create histograms and scatterplots

wavespeed1=sqrt(9.8*nccell1.bottomDepth.data);
wavespeed2=sqrt(9.8*nccell2.bottomDepth.data);

n=3
for col=1:3

ridgepack_multiplot(1,3,1,col)

if col==1
 x=wavespeed1; 
 y=wavespeed2; 
elseif col==2
 x=wavespeed1(latmasknorth); 
 y=wavespeed2(latmasknorth); 
elseif col==3
 x=wavespeed1(latmasksouth); 
 y=wavespeed2(latmasksouth); 
end

[xk,I] = sort(x);
yk = y(I);
p=polyfit(xk,yk,n)
yk=polyval(p,xk);

plot(x,y,'.','Color',0.85*[1 1 1])
h1=plot(xk,yk,'r')
h2=plot(xk,xk,'b')
xlabel('Unsmoothed (m/s)')
if col==1
 maxax=5+max([x;y])
 ylabel('100km-Smoothed (m/s)')
 title('Global Mesh','FontWeight','normal')
elseif col==2
 title('Polar: North of 60^{\circ}N','FontWeight','normal')
 legend([h1,h2],'Cubic Fit','1:1 Reference','Location','NorthWest')
 legend('boxoff')
elseif col==3
 title('Polar: South of 60^{\circ}S','FontWeight','normal')
end
xlim([0 maxax])
ylim([0 maxax])
axis square

end

ridgepack_multialign(gcf,'Impact of smoothing on barotropic wave speed: IcoswISC30E3')

cd(plotloc)
ridgepack_fprint('png',['Icos_Barotropic_Wave_Speed_Smoothed'],1,2)

return
  
