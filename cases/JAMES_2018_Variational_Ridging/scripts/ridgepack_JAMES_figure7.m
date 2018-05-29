% ridgepack_JAMES_figure7 - Generates Figure 7 in JAMES Variational Ridging paper
% 
% This script generates Figure 7 from:
%
% Roberts, A.F., E.C. Hunke, S.M. Kamal, W.H. Lipscomb, C. Horvat, W. Maslowski (2018),
% Variational Method for Sea Ice Ridging in Earth System Models, Part I: Theory, 
% submitted to J. Adv. Model Earth Sy.
%
% Andrew Roberts, Naval Postgraduate School, April 2018 (afrobert@nps.edu)

% version check
[v d] = version;
if str2num(d(end-3:end))<2018
 warning('This script designed for MATLAB 2018a or more recent version')
end

clear
close all

% parameter settings of scheme
hc=ridgepack_astroconstants;
rho=hc.rhoi.const;  % density of ice (kg/m^3)
rhos=hc.rhos.const; % density of snow (kg/m^3)
rhow=hc.rhow.const; % density of seawater (kg/m^3)

cols=lines(10); % colors

fontsize=12;

maxy=0;

for i=1:4

 if i==1
  hF=2.0;
  hFs=0.3;
  strain=-1/3;
  porosity=0;
  legendtext{i}='a) Imporous baseline ($\phi_R\!=\!0, \epsilon_{R_I}\!=\!-1/3$)';
  col=cols(1,:);
  pattern='-';
  linewidth=1.0;
 elseif i==2
  hF=2.0;
  hFs=0.3;
  strain=-1/2;
  legendtext{i}='b) Imporous, increased strain ($\phi_R\!=\!0, \epsilon_{R_I}\!=\!-1/2$)';
  porosity=0;
  col=cols(2,:);
  pattern=':';
  linewidth=1.0;
 elseif i==3
  hF=2.0;
  hFs=0.3;
  strain=-1/3;
  porosity=0.2;
  legendtext{i}='c) Porous ($\phi_R\!=\!0.2, \epsilon_{R_I}\!=\!-1/3$)';
  col=cols(4,:);
  pattern=':';
  linewidth=1.0;
 elseif i==4
  hF=2.0;
  hFs=0.3;
  strain=0;
  porosity=0.2;
  legendtext{i}='d) Porous, no strain ($\phi_R\!=\!0.2, \epsilon_{R_I}\!=\!0$)';
  col=cols(5,:);
  pattern=':';
  linewidth=1.0;
 elseif i==5
  hF=2.0;
  hFs=0.3;
  strain=0;
  porosity=0.12;
  legendtext{i}='$\phi_R\!=\!0.2, \epsilon_{R_I}\!=\!0$';
  col=cols(6,:);
  pattern=':';
  linewidth=0.75;
 end

 hFf=((rhow-rho)*hF+(rhow-rhos)*hFs)/rhow;

 if (hFs>hFf)
  error('Snow thickness too great')
 end

 hFd=(rho*hF+rhos*hFs)/rhow;

 hRs=hFs;
 hR=hF/(strain+1);

 hRf=((rhow-rho)*hR+(rhow-rhos)*hRs)/rhow;
 hRd=(rho*hR+rhos*hRs)/rhow;

 Hs=hFf+2*sqrt((hRd/(1-porosity)-hFd)*(hRf/(1-porosity)-hFf));
 Hk=(2*hRd/(1-porosity))-hFd;

 Lsd=(Hs-hFf);
 Lkd=(Hk-hFd);

 step(1)=hF;
 g(1)=0;

 step(2)=hF;
 g(2)=2/(2*Lkd);

 step(3)=hF+Lkd-Lsd;
 g(3)=2/(2*Lkd);

 step(4)=hF+Lkd-Lsd;
 g(4)=1/(2*Lkd);

 step(5)=hF+Lkd+Lsd;
 g(5)=1/(2*Lkd);

 step(6)=hF+Lkd+Lsd;
 g(6)=0;

 area=g(3)*(Lkd-Lsd)+g(5)*(2*Lsd);

 maxy=max(maxy,ceil(step(6)));

 set(gca,'FontSize',fontsize,'Box','on');
 hold on
 h(i)=plot(step,g,pattern,'Color',col,'LineWidth',linewidth);

end

legend(h,legendtext,'FontSize',fontsize,'Interpreter','Latex')
legend('boxoff')
xlim([hF-0.2 maxy])
ylim([0 1.1])
xlabel('$h$ (m)','Interpreter','Latex')
ylabel('Ridge ice thickness distribution $g_R(h,\phi_R)$','Interpreter','Latex')

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

