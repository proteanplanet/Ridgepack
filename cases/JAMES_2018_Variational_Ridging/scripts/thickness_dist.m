
clear
clf

rho=917.0; % density of ice (kg/m^3)
rhos=330.0; % density of snow (kg/m^3)
rhow=1026.0; % density of seawater (kg/m^3)

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

 step(6)=hF+Lkd+Lsd
 g(6)=0;

 area=g(3)*(Lkd-Lsd)+g(5)*(2*Lsd)

 maxy=max(maxy,ceil(step(6)));

 set(gca,'FontSize',fontsize,'Box','on');
 hold on
 h(i)=plot(step,g,pattern,'Color',col,'LineWidth',linewidth);

end

%set(gca,'YScale','log')

legend(h,legendtext,'FontSize',fontsize,'Interpreter','Latex')
legend('boxoff')
xlim([hF-0.2 maxy])
ylim([0 1.1])
xlabel('$h$ (m)','Interpreter','Latex')
ylabel('Ridge ice thickness distribution $g_R(h,\phi_R)$','Interpreter','Latex')

cd /Users/aroberts/Publications/2015_Unified_Morphology_1/figures

ncfprint('png',['Thickness_DistR'],1,2)
ncfprint('epsc',['Thickness_DistR'],1,2)
ncfprint('png',['Thickness_DistR_lowres'],1,1)
ncfprint('epsc',['Thickness_DistR_lowres'],1,1)











