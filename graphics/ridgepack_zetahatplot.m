

% generate underlying data switch
%generate=true;
generate=false;

% determine resolution of plot
%lowres=true;
lowres=false;

% directory for read/write of underlying data
cd ~/science/publications/2015_Variational_Principle/Data/Ridge_Energetics

% set latex as default interpreter
set(0,'DefaultTextInterpreter','Latex')

% generate or extract data
if generate

 [hfi,strainp,pep,phip,alphap,HK,HS,LK,LS]=paper_ridge_expected_plane(lowres);

 if lowres
  save expected_plane_lowres
 else
  save expected_plane
 end

else

 if lowres
  load expected_plane_lowres
 else
  load expected_plane
 end

end

% place data into a netcdf structure
nc.attributes.title='expected path';
[nc]=ncadd(nc,'hfi',hfi,'h_{fi}',{'hfi'},'m');
[nc]=ncadd(nc,'strain',strainp,'$\epsilon$',{'strain'},'');
[nc]=ncadd(nc,'pe',pep,'Potential Energy Density',{'hfi','strain'},'J m^{-2}');
[nc]=ncadd(nc,'pef',pep,'Filtered Potential Energy Density',{'hfi','strain'},'J m^{-2}');
[nc]=ncadd(nc,'phi',phip,'porosity',{'hfi','strain'},'');
[nc]=ncadd(nc,'alpha',alphap,'angle of repose',{'hfi','strain'},'degrees');
[nc]=ncadd(nc,'Hk',HK,'Keel Draft',{'hfi','strain'},'m');
[nc]=ncadd(nc,'Lk',LK,'Keel Width',{'hfi','strain'},'m');
[nc]=ncadd(nc,'Hs',HS,'Sail Height',{'hfi','strain'},'m');
[nc]=ncadd(nc,'Ls',LS,'Sail Width',{'hfi','strain'},'m');

% determine observed range and shade only every second cell outside the range
HSmax=5.24*sqrt(hfi); % based on (Tucker et al. 1984)
HSmaxfilt=ones(size(pep));
for i=1:length(hfi)
 idx=find(min(abs(HSmax(i)-HS(i,:)))==abs(HSmax(i)-HS(i,:)));
 HSstrain(i)=strainp(idx);
 HSmaxfilt(i,idx:end)=NaN;
end
nc.pef.data(1:2:end,1:2:end)=nc.pe.data(1:2:end,1:2:end).*HSmaxfilt(1:2:end,1:2:end);
nc.pef.data(2:2:end,2:2:end)=nc.pe.data(2:2:end,2:2:end).*HSmaxfilt(2:2:end,2:2:end);

% close any open figures
close all

% plot the data
ncimage(nc,'hfi','strain','pef',{},{},[1:1:11],10^-1,10^4,'vertical','parula')
set(gca,'Xscale','log','Ydir','reverse')
xmin=min(hfi);
xmax=max(hfi);
xlim([xmin xmax])
ymax=0;
ymin=min(strainp);
ylim([ymin ymax])
hold on
drawnow

% plot LK contours
drawlk=true;
%drawlk=false;

if drawlk

 cont=[0.1 1 10 100];
 conttitle=[1];
 linecol=1*[1 1 1];
 for i=1:length(cont)
  clear contlinex contliney contlinexl contlineyl;
  count=0;
  countl=0;
  for k=1:length(nc.hfi.data)
   idx=find(min(abs(nc.Lk.data(k,:)-cont(i)))==abs(nc.Lk.data(k,:)-cont(i)));
   if idx>1 & idx<length(nc.strain.data) & nc.hfi.data(k)>xmin     
    count=count+1;
    contlinex(count)=nc.hfi.data(k);
    contliney(count)=nc.strain.data(idx);
   end
   if idx>1 & idx<length(nc.strain.data) & ...
      nc.hfi.data(k)>0.001 & nc.strain.data(idx)>-0.80 
    countl=countl+1;
    contlinexl(countl)=nc.hfi.data(k);
    contlineyl(countl)=nc.strain.data(idx);
   end
  end
  if count>0
   line(contlinex,contliney,'Color',linecol,'LineStyle','-')
  end
  if countl>0
   ht=nctext(contlinexl,contlineyl,[num2str(cont(i),'%5.1f'),' m'],...
             10,linecol,'left','bottom','none');
   if cont(i)==conttitle
    hc=nctext(contlinexl,contlineyl,['Keel Width, $L_K$'],...
             10,linecol,'left','top','none');
   end
  end
 end
 drawnow

end


% plot HK contours
drawhk=true;
%drawhk=false;

if drawhk

 % draw keel draft
 hold on
 cont=[0.1 1 10 100];
 conttitle=[1];
 linecol=0*[1 1 1];
 for i=1:length(cont)
  clear contlinex contliney contlinexl contlineyl;
  count=0;
  countl=0;
  for k=1:length(nc.hfi.data)
   idx=find(min(abs(nc.Hk.data(k,:)-cont(i)))==abs(nc.Hk.data(k,:)-cont(i)));
   if idx>1 & idx<length(nc.strain.data) & nc.hfi.data(k)>xmin
    count=count+1;
    contlinex(count)=nc.hfi.data(k);
    contliney(count)=nc.strain.data(idx);
   end
   if idx>1 & idx<length(nc.strain.data) & ...
      nc.hfi.data(k)>xmin & nc.strain.data(idx)>-0.80
    countl=countl+1;
    contlinexl(countl)=nc.hfi.data(k);
    contlineyl(countl)=nc.strain.data(idx);
   end
  end
  if count>0
   line(contlinex,contliney,'Color',linecol)
  end
  if countl>0
   ht=nctext(contlinexl,contlineyl,[num2str(cont(i),'%5.1f'),' m'],...
            10,linecol,'left','bottom','none');
   if cont(i)==conttitle
    hc=nctext(contlinexl,contlineyl,['Keel Draft, $H_K$'],...
             10,linecol,'left','top','none');
   end
  end
 end
 drawnow

end

% plot Phi contours
drawphi=true;
%drawphi=false;

if drawphi

 % draw poroisity
 hold on
 conttitle=[0.35];
 cont=[0.05:0.05:0.35];
 %linecol=cols(7,:);
 %linecol=[1 0 0];
 linecol=[0.25 0 0.75];
 for i=1:length(cont)
  clear contlinex contliney;
  count=0;
  for k=1:length(nc.hfi.data)
   idx=find(min(abs(nc.phi.data(k,:)-cont(i)))==abs(nc.phi.data(k,:)-cont(i)));
   if idx>1 & idx<length(nc.strain.data) & nc.hfi.data(k)>xmin
    count=count+1;
    contlinex(count)=min(max(xmin,nc.hfi.data(k)),xmax);
    contliney(count)=nc.strain.data(idx);
   end
  end
  line([contlinex(1) contlinex(end)],...
      [contliney(1) contliney(end)],...
      'Color',linecol,'LineStyle','--')
  text(contlinex(end),contliney(end),[num2str(cont(i),'%5.2f'),'$\,\,\,\,\,$'],...
      'FontSize',10,'Color',linecol,...
      'HorizontalAlignment','right',...
      'VerticalAlignment','bottom','margin',1);
  if cont(i)==conttitle
   text(contlinex(end),contliney(end),['Porosity, $\phi_R \,\,\,\,\,$'],...
      'FontSize',10,'Color',linecol,...
      'HorizontalAlignment','right',...
      'VerticalAlignment','top','margin',1);
  end
 end
 drawnow

end

% plot paths
pathplot=true;
%pathplot=false;

if pathplot

 % path at 2m matching manifold diagram
 y=[0 -0.92];
 x=[2 2];
 plot(x,y,'r-')

 % arrow head of mean path
 text(x(end),y(end),...
                 ['$\bigtriangleup$'],...
                 'Color','r',...
                 'Fontsize',7,...
                 'Rotation',0,...
                 'VerticalAlignment','baseline',...
                 'HorizontalAlignment','center',...
                 'Interpreter','Latex')

 % overlay points from thickness distribution
 colplo={[0.5 0 1.0],[0.5 0 0.5],[0.75 0.25 0.25],[1 0.25 0]};
 erplor=[0,-0.2,-0.4,-0.6];
 for m=1:length(erplor)
   plot(x(1),erplor(m),'o','MarkerFaceColor','w','Color',colplo{m},'MarkerSize',4)
   plot(x(1),erplor(m),'+','Color',colplo{m},'MarkerSize',4)
 end

 % path at 0.5m matching manifold diagram
 x=[0.5 0.5];
 plot(x,y,'r:')

 % arrow head of mean path
 text(x(end),y(end),...
                 ['$\bigtriangleup$'],...
                 'Color','r',...
                 'Fontsize',7,...
                 'Rotation',0,...
                 'VerticalAlignment','baseline',...
                 'HorizontalAlignment','center',...
                 'Interpreter','Latex')

end

% label axes and title
ylabel('Strain, $\epsilon_R$')
xlabel('Parent Sheet Ice Thickness, $h_{F}$ (m)')
title('')
drawnow



% move to appropriate directory and write output
cd /Users/aroberts/Publications/2015_Unified_Morphology_1/figures
if lowres
 ncfprint('png',['Expected_Plane_lowres'],1,1)
 ncfprint('epsc',['Expected_Plane_lowres'],1,1)
else
 ncfprint('png',['Expected_Plane_midres'],1,1)
 ncfprint('epsc',['Expected_Plane_midres'],1,1)
 ncfprint('png',['Expected_Plane'],1,2)
 ncfprint('epsc',['Expected_Plane'],1,2)
end

