% ridgepack_JAMES_figure14 - Generates Figure 14 in JAMES Variational Ridging paper
% 
% This script generates Figure 14 from:
%
% Roberts, A.F., E.C. Hunke, S.M. Kamal, W.H. Lipscomb, C. Horvat, W. Maslowski (2018),
% Variational Method for Sea Ice Ridging in Earth System Models, Part I: Theory, 
% submitted to J. Adv. Model Earth Sy.
%
% This function generates and plots a zeta-hat plane as a test of the function 
% ridgepack_zetahatplane in the morphology library. 
%
% Andrew Roberts, Naval Postgraduate School, March 2018 (afrobert@nps.edu)

clear
close all

% generate or extract data on the zeta-hat plane (or not)
generate=true;
%generate=false;

% plot individual paths (or not)
pathplot=true;
%pathplot=false;

% resolution of phi and epsilon grid
%resolution=0.01; % lowest resolution 
%resolution=0.005; % moderate resolution
resolution=0.001; % high resolution

% set latex as default interpreter for graphics
set(0,'DefaultTextInterpreter','Latex')

% determine directory for read/write of zeta-hat plane data
dir=fileparts(which(mfilename));
writedir=[dir(1:strfind(dir,'scripts')-1),'output'];
[status,msg]=mkdir(writedir);
cd(writedir);

% generate or extract data for zeta-hat plane
if generate
 disp(['Generating data in: ',char(13),' ',writedir])
 [HF,EPSILON,PHI,ALPHAHAT,VR,HK,HS,LK,LS]=ridgepack_zetahatplane(resolution);
 save('ridgepack_zetahatplane','HF','EPSILON','PHI','ALPHAHAT','VR','HK','HS','LK','LS')
else
 disp(['Reading from:',char(13),' ',writedir])
 try
  load ridgepack_zetahatplane
 catch
  disp('Unable to load zeta-hat plane - switch to generate')
 end
end

% place data into a netcdf structure
nc.attributes.title='expected path';
[nc]=ridgepack_add(nc,'hf',HF,'h_{fi}',{'hf'},'m');
[nc]=ridgepack_add(nc,'epsilon',EPSILON,'$\epsilon$',{'epsilon'},'');
[nc]=ridgepack_add(nc,'phi',PHI,'porosity',{'hf','epsilon'},'');
[nc]=ridgepack_add(nc,'alpha',ALPHAHAT,'angle of repose',{'hf','epsilon'},'degrees');
[nc]=ridgepack_add(nc,'VR',VR,'Potential Energy Density',{'hf','epsilon'},'J m^{-2}');
[nc]=ridgepack_add(nc,'VRF',VR,'Filtered Potential Energy Density',{'hf','epsilon'},'J m^{-2}');
[nc]=ridgepack_add(nc,'Hk',HK,'Keel Draft',{'hf','epsilon'},'m');
[nc]=ridgepack_add(nc,'Lk',LK,'Keel Width',{'hf','epsilon'},'m');
[nc]=ridgepack_add(nc,'Hs',HS,'Sail Height',{'hf','epsilon'},'m');
[nc]=ridgepack_add(nc,'Ls',LS,'Sail Width',{'hf','epsilon'},'m');

% determine observed range and shade only every second cell outside the range
% based on Tucker et al. (1984)
HSmax=5.24*sqrt(HF); 
HSmaxfilt=ones(size(PHI));
for i=1:length(HF)
 idx=find(min(abs(HSmax(i)-HS(i,:)))==abs(HSmax(i)-HS(i,:)));
 HSepsilon(i)=EPSILON(idx);
 HSmaxfilt(i,idx:end)=NaN;
end
nc.VRF.data(1:2:end,1:2:end)=nc.VR.data(1:2:end,1:2:end).*HSmaxfilt(1:2:end,1:2:end);
nc.VRF.data(2:2:end,2:2:end)=nc.VR.data(2:2:end,2:2:end).*HSmaxfilt(2:2:end,2:2:end);

% close any open figures
figure(1)

% plot the data
ridgepack_image(nc,'hf','epsilon','VRF',{},{},[1:1:11],10^-1,10^4,'vertical','parula')
set(gca,'Xscale','log','Ydir','reverse')
xmin=min(HF);
xmax=max(HF);
xlim([xmin xmax])
ymax=0;
ymin=min(EPSILON);
ylim([ymin ymax])
hold on
drawnow

% plot LK contours
cont=[0.1 1 10 100];
conttitle=[1];
linecol=1*[1 1 1];
for i=1:length(cont)
  clear contlinex contliney contlinexl contlineyl;
  count=0;
  countl=0;
  for k=1:length(nc.hf.data)
   idx=find(min(abs(nc.Lk.data(k,:)-cont(i)))==abs(nc.Lk.data(k,:)-cont(i)));
   if idx>1 & idx<length(nc.epsilon.data) & nc.hf.data(k)>xmin     
    count=count+1;
    contlinex(count)=nc.hf.data(k);
    contliney(count)=nc.epsilon.data(idx);
   end
   if idx>1 & idx<length(nc.epsilon.data) & ...
      nc.hf.data(k)>0.001 & nc.epsilon.data(idx)>-0.80 
    countl=countl+1;
    contlinexl(countl)=nc.hf.data(k);
    contlineyl(countl)=nc.epsilon.data(idx);
   end
  end
  if count>0
   line(contlinex,contliney,'Color',linecol,'LineStyle','-')
  end
  if countl>0
   ht=ridgepack_text(contlinexl,contlineyl,[num2str(cont(i),'%5.1f'),' m'],...
             10,linecol,'left','bottom','none');
   if cont(i)==conttitle
    hc=ridgepack_text(contlinexl,contlineyl,['Keel Width, $L_K$'],...
             10,linecol,'left','top','none');
   end
  end
end
drawnow

% plot HK contours
hold on
cont=[0.1 1 10 100];
conttitle=[1];
linecol=0*[1 1 1];
for i=1:length(cont)
  clear contlinex contliney contlinexl contlineyl;
  count=0;
  countl=0;
  for k=1:length(nc.hf.data)
   idx=find(min(abs(nc.Hk.data(k,:)-cont(i)))==abs(nc.Hk.data(k,:)-cont(i)));
   if idx>1 & idx<length(nc.epsilon.data) & nc.hf.data(k)>xmin
    count=count+1;
    contlinex(count)=nc.hf.data(k);
    contliney(count)=nc.epsilon.data(idx);
   end
   if idx>1 & idx<length(nc.epsilon.data) & ...
      nc.hf.data(k)>xmin & nc.epsilon.data(idx)>-0.80
    countl=countl+1;
    contlinexl(countl)=nc.hf.data(k);
    contlineyl(countl)=nc.epsilon.data(idx);
   end
  end
  if count>0
   line(contlinex,contliney,'Color',linecol)
  end
  if countl>0
   ht=ridgepack_text(contlinexl,contlineyl,[num2str(cont(i),'%5.1f'),' m'],...
            10,linecol,'left','bottom','none');
   if cont(i)==conttitle
    hc=ridgepack_text(contlinexl,contlineyl,['Keel Draft, $H_K$'],...
             10,linecol,'left','top','none');
   end
  end
end
drawnow

% plot Phi contours
hold on
conttitle=[0.35];
cont=[0.05:0.05:0.35];
linecol=[0.25 0 0.75];
for i=1:length(cont)
  clear contlinex contliney;
  count=0;
  for k=1:length(nc.hf.data)
   idx=find(min(abs(nc.phi.data(k,:)-cont(i)))==abs(nc.phi.data(k,:)-cont(i)));
   if idx>1 & idx<length(nc.epsilon.data) & nc.hf.data(k)>xmin
    count=count+1;
    contlinex(count)=min(max(xmin,nc.hf.data(k)),xmax);
    contliney(count)=nc.epsilon.data(idx);
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

% plot individual paths paths
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
ylabel('Strain, $\epsilon_{R_I}$')
xlabel('Parent Sheet Ice Thickness, $h_{F}$ (m)')
title('')

% determine directory for read/write
dir=fileparts(which(mfilename));
outdir=[dir(1:strfind(dir,'scripts')-1),'output'];
[status,msg]=mkdir(outdir);
cd(outdir);

% determine filename
x=strfind(mfilename,'_');
thisfilename=mfilename;
graphicsout=thisfilename(x(end)+1:end);

% output
disp(['Writing graphics output ',graphicsout,' to:',char(13),' ',pwd])

% print figure
ridgepack_fprint('epsc',graphicsout,1,2)
ridgepack_fprint('png',graphicsout,1,2)

