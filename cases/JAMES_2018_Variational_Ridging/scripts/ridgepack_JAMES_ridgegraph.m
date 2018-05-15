% ridgepack_JAMES_ridgegraph - Generates ridgegraphs JAMES Variation Ridging paper
% 
% This script generates ridgegraphs from:
%
% Roberts, A.F., E.C. Hunke, S.M. Kamal, W.H. Lipscomb, C. Horvat, W. Maslowski (2018),
% Variational Method for Sea Ice Ridging in Earth System Models, Part I: Theory, 
% submitted to J. Adv. Model Earth Sy.
%
% The following settings generate figures in the paper:
%
% fig=1 - Generates figure 8 
% fig=2 - Generates figure 9 
% fig=3 - Generates figure S1 
% fig=6 - Generates figure S2 
%
% Andrew Roberts, Naval Postgraduate School, April 2018 (afrobert@nps.edu)

clear
close all

fig=1

if fig==1
 nrows=1;
 ncols=4;
elseif fig==2
 nrows=1;
 ncols=2;
elseif fig==3
 nrows=3;
 ncols=1;
elseif fig==4
 nrows=3;
 ncols=1;
elseif fig==5
 nrows=3;
 ncols=1;
elseif fig==6
 nrows=1;
 ncols=2;
elseif fig==7
 nrows=1;
 ncols=1;
end

% Create figure
h=figure(1);
%set(h,'Position',[0 0 ncols*250 nrows*250/2]);
clf

% geometric settings of figure
if fig==1
 boxleftx=-0.06
elseif fig==2
 boxleftx=+0.13
elseif fig==3
 boxleftx=+0.10;
elseif fig==4
 boxleftx=+0.14;
elseif fig==5
 boxleftx=+0.17;
elseif fig==6
 boxleftx=+0.15;
elseif fig==7
 boxleftx=+0.15;
else 
 boxleftx=0.05;
end
if fig==6
 boxrightx=0.72;
else
 boxrightx=1-boxleftx;
end
boxcenterx=(boxleftx+boxrightx)/2;
figLk=boxrightx-boxleftx;
levelextent=0.0;
if fig==1
 sealeft=-0.045
 searight=1.255
elseif fig==2
 sealeft=0.03
 searight=0.91
else
 sealeft=0.0
 searight=1.0
end
aspectratio=1;
textoffset=0.012;
arrowhead=0.010;
envelope=0.025
bluecol=[0 0.447058826684952 0.74117648601532];
greycol=0.25*[1 1 1];

% parameter settings of scheme
hc=ridgepack_astroconstants;
rho=hc.rhoi.const;  % density of ice (kg/m^3)
rhos=hc.rhos.const; % density of snow (kg/m^3)
rhow=hc.rhow.const; % density of seawater (kg/m^3)
g=hc.ghat.const; % density of seawater (kg/m^3)

if fig==1 
 notation='bacd';
elseif fig==3 | fig==4 | fig==5
 notation='bac';
else
 notation='abcdefgh';
end

if fig==4
 cd /Users/aroberts/science/publications/2015_Variational_Principle/Data/Ridge_Energetics
 load meanpath1
end

if nrows==1 | ncols==1
 maxset=max(ncols,nrows)
else
 error('not set up for multiple rows AND multiple columns')
end

for setting=1:maxset

 % level ice and snow thickness
 hli=2.0; % level ice thickness
 hls=0.3; % snow thickness

 % deformed ice and snow thickness
 hdi=3.0; % ridged ice thickness
 hds=0.3; % ridged snow thickness
 
 alphar=22; % angle of ridge
 alphak=22; % angle of keel

 porosity=0.2; % porosity of ridge and keel complex

 cols=lines(10); % set line colors

 if fig==1
  if setting==1
   edgecol=cols(2,:);
   epsilon=-1/2;
   hdi=hli/(epsilon+1);
   porosity=0.; % porosity of ridge and keel complex
   alpha_R=22; %angle of compressional repose
   theta=0; %angle of actual repose
  elseif setting==2
   edgecol=cols(1,:);
   epsilon=-1/3;
   hdi=hli/(epsilon+1);
   porosity=0; % porosity of ridge and keel complex
   alpha_R=22; %angle of compressional repose
   theta=0; %angle of actual repose
  elseif setting==3
   edgecol=cols(4,:);
   epsilon=-1/3;
   hdi=hli/(epsilon+1);
   porosity=0.2; % porosity of ridge and keel complex
   alpha_R=22; %angle of compressional repose
   theta=0; %angle of actual repose
  elseif setting==4
   edgecol=cols(5,:);
   epsilon=0;
   hdi=hli/(epsilon+1);
   porosity=0.2; % porosity of ridge and keel complex
   alpha_R=22; %angle of compressional repose
   theta=0; %angle of actual repose
  end
  alphar=acotd(cosd(theta)*cotd(alpha_R)); % apparent angle of repose
  alphak=alphar % angle of keel
 elseif fig==2
  if setting==1
   edgecol=cols(4,:);
   epsilon=-1/3;
   hdi=hli/(epsilon+1);
   porosity=0.2; % porosity of ridge and keel complex
   alpha_R=22; %angle of compressional repose
   theta=0; %angle of actual repose
  elseif setting==2
   edgecol=cols(4,:);
   epsilon=-1/3;
   hdi=hli/(epsilon+1);
   porosity=0.2; % porosity of ridge and keel complex
   alpha_R=22; %angle of compressional repose
   theta=55; %angle of actual repose
  end
  alphar=acotd(cosd(theta)*cotd(alpha_R)); % apparent angle of repose
  alphak=alphar % angle of keel
 elseif fig==3
  edgecol=[0 0 0];
  if setting==1
   hli=2.0
   hls=0.3
   hdi=3.0
   epsilon=(hli-hdi)/hdi
   hds=0.5*hls/(epsilon+1)
  elseif setting==2
   hls=0.0
   hli=2.0
   hds=0.0
   hdi=3.0
  elseif setting==3
   hls=0.3
   hli=2.0
   hds=0.3
   hdi=3.0
  end
 elseif fig==4
  if setting==1
   edgecol='g';
   hls=hfsP1;
   hli=hfiP1;
   hds=hfsP1;
   hdi=hfiP1./(1-stiP1);
   alphak=alphaP1;
   alphar=alphaP1;
   porosity=phiP1;
  elseif setting==2
   edgecol='b';
   hls=hfsP2;
   hli=hfiP2;
   hds=hfsP2;
   hdi=hfiP2./(1-stiP2);
   alphak=alphaP2;
   alphar=alphaP2;
   porosity=phiP2;
  elseif setting==3
   edgecol='m';
   hls=hfsP3;
   hli=hfiP3;
   hds=hfsP3;
   hdi=hfiP3./(1-stiP3);
   alphak=alphaP3;
   alphar=alphaP3;
   porosity=phiP3;
  end
 elseif fig==5
  edgecol=[0 0 0];
  if setting==1
   hls=0.3
   hli=2.0
   hds=0.3
   hdi=3.0
   alphak=26.6;
   alphar=26.6;
  elseif setting==2
   hls=0.3
   hli=2.0
   hds=0.3
   hdi=3.0
   alphak=26.6;
   alphar=32.9;
  elseif setting==3
   hls=0.3
   hli=2.0
   hds=0.3
   hdi=3.0
   alphak=26.6;
   alphar=20.7;
  end
 elseif fig==6
  edgecol=[0 0 0];
  if setting==1
   hls=0.3
   hli=2.0
   hds=0.3
   hdi=3.0
   alphak=22;
   alphar=12;
  elseif setting==2
   hls=0.3
   hli=2.0
   hds=0.3
   hdi=3.0
   alphak=22;
   alphar=32;
  end
 elseif fig==7
  edgecol=[0 0 0];
  epsilon=-1/3;
  hdi=hli/(epsilon+1);
  porosity=0; % porosity of ridge and keel complex
  alpha_R=22; %angle of compressional repose
  theta=0; %angle of actual repose
 else
  error('fig is wrong')
 end

 % calculate freeboard and draft of level ice
 hld=(rho*hli+rhos*hls)/rhow; % level draft
 hlf=(hli+hls)-hld; % level freeboard

 % calculate freeboard and draft of deformed ice 
 hdd=(rho*hdi+rhos*hds)/rhow; % ridged draft
 hdf=(hdi+hds)-hdd; % ridged freeboard

 if hlf<hls
  error('Level snow too thick')
 elseif hdf<hds
  error('Deformed snow too thick')
 end

 % check for bounds of level ice
 if hld<0
  hld=0;
  hlf=0;
 elseif hld>hdd/(1-porosity) | hlf>hdf/(1-porosity)
  hdd=hld;
  hdf=hlf;
 end

 % calculate depth of keel relative to sea level
 Hk=(2*hdd/(1-porosity))-hld;

 % calculate horizontal extent of keel structure 
 Lk=2*(Hk-hld)/tan(alphak*pi/180)

 % calculate height of ridge
 Hr=hlf+sqrt((tan(alphar*pi/180))*(((hdf*Lk)/(1-porosity))-hlf*Lk));

 % calculate horizontal extent of ridge structure 
 Lr=2*(Hr-hlf)/tan(alphar*pi/180)

 % use this to determine the scale factor for the diagram
 if setting==1; 
  scalefactorx=figLk/Lk;
  scalefactory=min(aspectratio*scalefactorx,1/(Hk+Hr));
  Lk1=Lk;
  Hk1=Hk;
  Hr1=Hr;
 end

 if fig==5 % calculate form drag
  Hrd=Hr-hlf;
  Scsquared=(1-exp(-0.18*Lk/Hrd));
  Cra=0.2;
  Conc=1;
  z0=5*(10^-4);
  Cran=0.5*Cra*Scsquared*(Hrd/Lk)*Conc*(log(Hrd/z0)/log(20/z0))^2
 end

 % calculate level ice box and surface snow, centered slightly to left
 sealevely=0.6; % rho/rhow
 ylsnow=sealevely+hlf*scalefactory
 ylbottom=sealevely-hld*scalefactory;
 xlleft=boxleftx+(boxrightx-boxleftx-Lk*scalefactory)/2
 xlright=boxrightx-(boxrightx-boxleftx-Lk*scalefactory)/2

 % calculate ridge top and keel bottom
 yrtop=ylsnow+(Hr-hlf)*scalefactory;
 ykbottom=ylbottom-(Hk-hld)*scalefactory;

 if fig==8
  yrtop1=yrtop;
  ykbottom1=ykbottom;
  xlright1=xlright;
 elseif setting==1
  yrtop1=yrtop;
  ykbottom1=ykbottom;
  xlright1=xlright;
 end

 if fig==1 | fig==2
  xlright1=xlright;
 end

 % set axis location order by setting
 if fig==1 
  row(1)=1; col(1)=2;
  row(2)=1; col(2)=1;
  row(3)=1; col(3)=3;
  row(4)=1; col(4)=4;
  if setting==1
   offset=-0.00;
  elseif setting==2
   offset=-0.05;
  elseif setting==3
   offset=+0.07;
  elseif setting==4
   offset=+0.04;
  end
 elseif fig==2 
  row(1)=1; col(1)=1;
  row(2)=1; col(2)=2;
  if setting==1
   offset=0.05;
  elseif setting==2
   offset=-0.01;
  end
 elseif fig==3
  row(1)=2; col(1)=1;
  row(2)=1; col(2)=1;
  row(3)=3; col(3)=1;
  offset=0.0;
 elseif fig==4
  row(1)=1; col(1)=1;
  row(2)=2; col(2)=1;
  row(3)=3; col(3)=1;
  offset=0.0;
 elseif fig==5
  row(1)=1; col(1)=1;
  row(2)=2; col(2)=1;
  row(3)=3; col(3)=1;
  offset=0.0;
 elseif fig==6
  row(1)=1; col(1)=1;
  row(2)=1; col(2)=2;
  offset=0.0;
 elseif fig==7
  row(1)=1; col(1)=1;
  offset=0.0;
 end

 % start axis
 axes('Position',[(col(setting)-1)/ncols+offset ...
                  (nrows-row(setting))/nrows 1/ncols 1/nrows],...
      'DataAspectRatio',[1 1 1],'PlotBoxAspectRatio',[1 1 1],...
      'Xlim',[0 1],'Ylim',[0.25 0.75],'ActivePositionProperty','position',...
      'Box','off','Visible','off','Clipping','off','Xtick',[],'XTickLabel',[],...
      'Ytick',[],'YTickLabel',[])
 
 ridgepack_clearax

 if fig~=7
  % add "a)", "b)", and "c)" etc
  ymin=get(gca,'Ylim');
  if nrows>1 & ncols==1
   text(0,sealevely+envelope,[notation(setting),')'],...
          'FontName','Helvetica',...
          'Interpreter','Tex')
  else
   text(xlleft-envelope,sealevely+Hr1*scalefactory,[notation(setting),')'],...
          'FontName','Helvetica',...
          'Interpreter','Tex')
  end
 end

 set(0,'DefaultTextInterpreter','Latex')

 if fig==2 | fig==6 | fig==3
  set(gca,'FontSize',12)
 elseif fig==1
  set(gca,'FontSize',9)
 else
  set(gca,'FontSize',11)
 end
 labelsize=get(gca,'FontSize');

 % fill in ridge and freeboard
 x=[xlleft, xlleft, boxcenterx-scalefactorx*Lr/2, boxcenterx,...
     boxcenterx+scalefactorx*Lr/2, xlright, xlright, boxcenterx, xlleft];
 y=[ylbottom, ylsnow, ylsnow, yrtop, ylsnow, ylsnow, ylbottom, ykbottom, ylbottom];
 if fig==7
  patch(x,y,0.95*[1 1 1],'EdgeColor',edgecol,'FaceColor','none')
 else
  patch(x,y,0.95*[1 1 1],'EdgeColor',edgecol)
 end

 % plot sea level line across figure
 x=[sealeft searight]; 
 y=[sealevely sealevely];
 line(x,y,'Color',bluecol,'LineStyle','-.');
 if fig==7
  text(boxcenterx,sealevely,...
           {'sea level'},...
           'VerticalAlignment','bottom',...
           'HorizontalAlignment','center',...
           'Fontsize',labelsize,...
           'Color',bluecol,...
           'EdgeColor','none')
 elseif (setting==1 & fig~=3) | (setting==2 & fig==3) 
  text(boxcenterx,sealevely-0.75*textoffset,...
           {'sea level'},...
           'VerticalAlignment','bottom',...
           'HorizontalAlignment','center',...
           'Fontsize',labelsize-0.25,...
           'Color',bluecol,...
           'EdgeColor','none')
 end

 % plot relevant metric information
 if fig==1
   if epsilon<0
    epstring=num2str(epsilon,'%8.2f');
   else
    epstring=num2str(epsilon);
   end
   if porosity>0
    porstring=num2str(porosity,'%8.2f');
   else
    porstring=num2str(porosity);
   end
   textbox={['\makebox[4in][c]{$\phi_R\!=\!',porstring,'$}'],...
            ['\makebox[4in][c]{$\epsilon_{R_I}\!=\!',epstring,'$}']}
 elseif fig==2
   textbox={['\makebox[4in][c]{$\alpha_R=',num2str(alpha_R,'%8.1f'),'^{\circ}$}'],...
            ['\makebox[4in][c]{$\theta_R=',num2str(180-theta,'%8.1f'),'^{\circ}$}'],...
            ['\makebox[4in][c]{$\alpha_K=',num2str(alphak,'%8.1f'),'^{\circ}$}']}
 elseif fig==3
   textbox={['\makebox[4in][c]{$h_{Fs}=',num2str(hls),'$ m}'],...
            ['\makebox[4in][c]{$h_{Rs}=',num2str(hds),'$ m}']}
 elseif fig==4
   textbox={['\makebox[4in][c]{$h_{Fi}=',num2str(hli),'$ m}'],...
            ['\makebox[4in][c]{$h_{Fs}=',num2str(hls),'$ m}'],...
            ['\makebox[4in][c]{$h_{Ri}=',num2str(hdi,'%5.2f'),'$ m}'],...
            ['\makebox[4in][c]{$\phi_R=',num2str(porosity),'$}']}
 elseif fig==5
   textbox={['\makebox[4in][c]{$\alpha_S=',num2str(alphar),'^{\circ}$}'],...
            ['\makebox[4in][c]{$L_S=',num2str(Lr,'%5.1f'),'$ m}'],...
            ['\makebox[4in][c]{$C_{Ran}=',num2str(Cran*(10^3),'%5.1f'),'\times 10^{-3}$}']}
 elseif fig==6
   textbox={['\makebox[4in][c]{$\alpha_S=',num2str(alphar),'^{\circ}$}'],...
            ['\makebox[4in][c]{$L_S=',num2str(Lr,'%5.1f'),'$ m}']}
 elseif fig==7
   textbox='';
 end

 % add text into this 
 text(boxcenterx,sealevely-1.5*textoffset,textbox,...
           'VerticalAlignment','cap','HorizontalAlignment','center',...
           'Fontsize',labelsize-0.25,'Color',edgecol,'EdgeColor','none')

 % add path indicators
 if fig==4
  if setting==1
   textbox='$\xi_1$';
   pcols='g';
  elseif setting==2
   textbox='$\xi_2$';
   pcols='b';
  elseif setting==3
   textbox='$\xi_3$';
   pcols='m';
  end
  text(xlleft+1.0*textoffset,sealevely-0.1*textoffset,textbox,...
           'VerticalAlignment','top',...
           'HorizontalAlignment','left',...
           'Fontsize',labelsize,...
           'Color',pcols,...
           'EdgeColor','none')
 end

 % provide horizontal scale bar
 x=[xlleft xlright];
 if fig==4
  y=[ykbottom-envelope ykbottom-envelope];
 else
  y=[ykbottom1-envelope ykbottom1-envelope];
 end
 if fig==7 
  textbox=['$L_{K}$'];
 elseif fig==3 
  textbox=['$L_{K}\!=\!',num2str(Lk,'%8.1f'),'$\,m'];
 elseif fig==1 & setting==2
  textbox=['$L_{K}\!=\!',num2str(Lk,'%8.1f'),'$\,m'];
 elseif setting==1
  textbox=['$L_{K}\!=\!',num2str(Lk,'%8.1f'),'$\,m'];
 else
  textbox=['$',num2str(Lk,'%8.1f'),'$\,m'];
 end
 line(x,y,'Color',0.0*[1 1 1],'LineStyle','-');
 line([x(1) x(1)],[y(1)-arrowhead y(1)+arrowhead],...
           'Color',0.0*[1 1 1],'LineStyle','-');
 line([x(2) x(2)],[y(1)-arrowhead y(1)+arrowhead],...
            'Color',0.0*[1 1 1],'LineStyle','-');
 text(sum(x)/2,y(1)-textoffset,textbox,...
           'VerticalAlignment','top',...
           'Fontsize',labelsize,...
           'HorizontalAlignment','center',...
           'EdgeColor','none')

 % provide ridge height scale bar
 x=[xlright1+envelope xlright1+envelope];
 y=[sealevely yrtop];
 line(x,y,'Color',0.0*[1 1 1],'LineStyle','-');
 handle=line([x(1)-arrowhead x(1)+arrowhead],[y(1) y(1)],...
          'Color',0.0*[1 1 1],'LineStyle','-');
 uistack(handle,'top')
 handle=line([x(1)-arrowhead x(1)+arrowhead],[y(2) y(2)],...
          'Color',0.0*[1 1 1],'LineStyle','-');
 uistack(handle,'top')
 if fig==1 & setting~=2 | (fig==2 | fig==6) & setting~=1 
  text(x(1)+textoffset,sum(y)/2,...
           [num2str(Hr,'%8.1f'),'\,m'],...
           'VerticalAlignment','middle',...
	   'Interpreter','Latex',...
           'Fontsize',labelsize,...
           'HorizontalAlignment','left',...
           'EdgeColor','none')
 elseif fig==7
  text(x(1)+textoffset,sum(y)/2,...
           ['$H_S$'],...
           'VerticalAlignment','middle',...
	   'Interpreter','Latex',...
           'Fontsize',labelsize,...
           'HorizontalAlignment','left',...
           'EdgeColor','none')
 
 else
  text(x(1)+textoffset,sum(y)/2,...
           ['$H_S\!=$',num2str(Hr,'%8.1f'),'\,m'],...
           'VerticalAlignment','middle',...
	   'Interpreter','Latex',...
           'Fontsize',labelsize,...
           'HorizontalAlignment','left',...
           'EdgeColor','none')
 end

 % provide keel depth scale bar
 x=[xlright1+envelope xlright1+envelope];
 y=[sealevely ykbottom];
 line(x,y,'Color',0.0*[1 1 1],'LineStyle','-');
 handle=line([x(1)-arrowhead x(1)+arrowhead],[y(1) y(1)],...
          'Color',0.0*[1 1 1],'LineStyle','-');
 uistack(handle,'top')
 handle=line([x(1)-arrowhead x(1)+arrowhead],[y(2) y(2)],...
          'Color',0.0*[1 1 1],'LineStyle','-');
 uistack(handle,'top')
 if fig==1 & setting~=2 | (fig==2 | fig==6) & setting~=1 
  text(x(1)+textoffset,sum(y)/2,...
           [num2str(Hk,'%8.1f'),'\,m'],...
           'VerticalAlignment','middle',...
	   'Interpreter','Latex',...
           'Fontsize',labelsize,...
           'HorizontalAlignment','left',...
           'EdgeColor','none')
 elseif fig==7
  text(x(1)+textoffset,sum(y)/2,...
           ['$H_K$'],...
           'VerticalAlignment','middle',...
	   'Interpreter','Latex',...
           'Fontsize',labelsize,...
           'HorizontalAlignment','left',...
           'EdgeColor','none')
 else
  text(x(1)+textoffset,sum(y)/2,...
           ['$H_{K}\!=$',num2str(Hk,'%8.1f'),'\,m'],...
           'VerticalAlignment','middle',...
	   'Interpreter','Latex',...
           'Fontsize',labelsize,...
           'HorizontalAlignment','left',...
           'EdgeColor','none')
 end


end % for setting=1:3

% determine directory for read/write
dir=fileparts(which(mfilename));
cd([dir(1:strfind(dir,'scripts')-1),'output']);

% determine filename
if fig==1
 graphicsout='figure8';
elseif fig==2
 graphicsout='figure9';
elseif fig==3
 graphicsout='figureS1';
elseif fig==6
 graphicsout='figureS2';
else
 graphicsout=['Comparison_Diagram_',num2str(fig)];
end

% output
disp(['Writing graphics output ',graphicsout,' to ',pwd])

% print figure
ridgepack_fprint('epsc',graphicsout,1,2)


