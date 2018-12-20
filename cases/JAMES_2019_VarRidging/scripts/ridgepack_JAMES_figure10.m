% ridgepack_JAMES_figure10 - Generates Figure 10 in JAMES Variational Ridging paper
% 
% This script generates Figure 10 from:
%
% Roberts, A.F., E.C. Hunke, S.M. Kamal, W.H. Lipscomb, C. Horvat, 
% W. Maslowski (2019), Variational Method for Sea Ice Ridging in 
% Earth System Models, J. Adv. Model Earth Sy.
%
% VERSION/CASES: Ridgepack 1.0.1/JAMES_2019_VarRidging
%
% CONTACT: Andrew Roberts, afroberts@lanl.gov 
%
% FILE HISTORY:
% Author: Andrew Roberts, Naval Postgraduate School, April 2018 
% Update: Andrew Roberts, Los Alamos National Laboratory, December 2018

% version check
[v d] = version;
if str2num(d(end-3:end))<2018
 warning('This script designed for MATLAB 2018a or more recent version')
end

clear
close all

% Create figure
h=figure(1);

% set up colors
cols=lines(10);
cols(11,:)=[0 0 0];
cols(11,:)=[0 0 1];

boxleftx=0.002;
boxrightx=1-boxleftx;

figLk=boxrightx-boxleftx;
levelextent=0.0;
sealeft=-0.05;
searight=1.05;
aspectratio=1;
textoffset=0.012;
arrowhead=0.005;
envelope=0.025;
bluecol=[0 0.447058826684952 0.74117648601532];
gcol=cols(7,:);
greycol=0.0*[1 1 1];
icecol=0.75*[1 1 1];
bluecol=[0 0.447058826684952 0.74117648601532];
whitecol=0.999999*[1 1 1];
fricol=11;

% parameter settings of scheme
hc=ridgepack_astroconstants;
rho=hc.rhoi.const;  % density of ice (kg/m^3)
rhos=hc.rhos.const; % density of snow (kg/m^3)
rhow=hc.rhow.const; % density of seawater (kg/m^3)

% level ice and snow thickness left side
hfi1=2.0; % level ice thickness
hfs1=0.0; % snow thickness

% level ice and snow thickness right side
hfi2=1.0; % level ice thickness
hfs2=0.0; % snow thickness

alpha1=18; % angle of ridge/keel repose
alpha2=24; % angle of ridge/keel repose

% deformed ice and snow thickness
hdi1=3.0; % ridged ice thickness
hds1=0.0; % ridged snow thickness

% porosity of ridge and keel complex
porosity=0.2; 

% calculate freeboard and draft of level ice
hfd1=(rho*hfi1+rhos*hfs1)/rhow; % level draft
hff1=(hfi1+hfs1)-hfd1; % level freeboard

% calculate freeboard and draft of level ice
hfd2=(rho*hfi2+rhos*hfs2)/rhow; % level draft
hff2=(hfi2+hfs2)-hfd2; % level freeboard

% calculate freeboard and draft of deformed ice 
hdd1=(rho*hdi1+rhos*hds1)/rhow; % ridged draft
hdf1=(hdi1+hds1)-hdd1; % ridged freeboard

% check for bounds of level ice
if hfd1<0
 hfd1=0;
 hff1=0;
elseif hfd1>hdd1/(1-porosity) | hff1>hdf1/(1-porosity)
 hdd1=hfd1;
 hdf1=hff1;
end

% calculate depth of keel relative to sea level
Hk=(2*hdd1/(1-porosity))-hfd1;

% calculate horizontal extent of keel structure 
Lk1=2*(Hk-hfd1)/tan(alpha1*pi/180);

% calculate horizontal extent of keel structure 
Lk2=2*(Hk-hfd2)/tan(alpha2*pi/180);

% calculate height of ridge
Hr=hff1+sqrt((tan(alpha1*pi/180))*(((hdf1*Lk1)/(1-porosity))-hff1*Lk1));

% calculate horizontal extent of ridge structure 
Lr1=2*(Hr-hff1)/tan(alpha1*pi/180);

% calculate horizontal extent of ridge structure 
Lr2=2*(Hr-hff2)/tan(alpha2*pi/180);

% position of sea level reference
sealevely=0.6; 

% use this to determine the scale factor for the diagram
scalefactorx=figLk/Lk1;
scalefactory=min(aspectratio*scalefactorx,1/(Hk+Hr));

% locate box center
boxcenterxt=(boxleftx*Lk2+boxrightx*Lk1)./(Lk1+Lk2);
boxcenterxb=(boxleftx*Lk2+boxrightx*Lk1)./(Lk1+Lk2);

% calculate level ice box and surface snow, centered slightly to left
ylsnow1=sealevely+hff1*scalefactory;
ylbottom1=sealevely-hfd1*scalefactory;
xlleft=(boxleftx*Lk2+boxrightx*Lk1)./(Lk1+Lk2)-(Lk1*scalefactorx)/2;

ylsnow2=sealevely+hff2*scalefactory;
ylbottom2=sealevely-hfd2*scalefactory;
xlright=(boxleftx*Lk2+boxrightx*Lk1)./(Lk1+Lk2)+(Lk2*scalefactorx)/2;

% calculate ridge top and keel bottom
yrtop=ylsnow1+(Hr-hff1)*scalefactory;
ykbottom=ylbottom1-(Hk-hfd1)*scalefactory;

% locate centroid
yb=sealevely-scalefactory*hdd1/(1-porosity);
yt=sealevely+scalefactory*(hdf1-hds1)/(1-porosity);
ys=sealevely+scalefactory*hdf1/(1-porosity);
yc=(0.5*(yt+yb)*rho + 0.5*(ys+yt)*rhos)./(rho+rhos);
xc=mean([xlleft xlright]);

% start axis
axes('Position',[0.025 0.025 0.95 0.95],...
     'DataAspectRatio',[1 1 1],'PlotBoxAspectRatio',[1 1 1],...
     'Xlim',[0 1],'Ylim',[0 1],'ActivePositionProperty','position',...
     'Box','off','Visible','off','Clipping','off','Xtick',[],'XTickLabel',[],...
     'Ytick',[],'YTickLabel',[])
 
ridgepack_clearax

ymin=get(gca,'Ylim');

set(gca,'FontSize',12)
labelsize=get(gca,'FontSize');

% calculate area above sea level in plotting coordinates for idealized shape
freeboardarea1=0.5*scalefactorx*Lr1*(yrtop-ylsnow1);
draftarea1=0.5*scalefactorx*Lk1*(ylbottom1-ykbottom);

freeboardarea2=0.5*scalefactorx*Lr2*(yrtop-ylsnow2);
draftarea2=0.5*scalefactorx*Lk2*(ylbottom2-ykbottom);

% plot roughened ridge structure with these parameters
maxx=25; % grid points along half-keel
jx=[2:5:24]; % randomly grab index to adjust
delt=3;  % size of delta to draw keels
if max(jx)>maxx
 error('maxx is too small')
end

% calculate left half keel
xincrement=(boxcenterxb-xlleft)./maxx;
yincrement=(ylbottom1-ykbottom)./maxx;
yincrementfirst=yincrement;

x=zeros([1 maxx+1]);
y=zeros([1 maxx+1]);

while any(y(2:end)>=ylbottom1 | y(1:end-1)<=ykbottom)
 area=zeros([maxx 1]);
 for i=1:maxx+1
  x(i)=xlleft+(i-1)*xincrement;
  coeff=(i/(maxx+1))^4;
  if i==1
   y(i)=ylbottom1;
  elseif i==maxx+1
   y(i)=ykbottom;
   area(i)=(x(i)-x(i-1))*(ylbottom1-0.5*(y(i)+y(i-1)));
  else
   y(i)=y(i-1)-delt*yincrement*(rand-0.49);
   y(i)=max(ykbottom+3*yincrement,min(y(i),ylbottom1-3*yincrement));
   y(i)=(1-coeff)*y(i)+coeff*ykbottom;
   area(i)=(x(i)-x(i-1))*(ylbottom1-0.5*(y(i)+y(i-1)));
  end
 end
 totalareabelow=sum(area(:));

 for idx=jx
  newarea=sum(area(idx:idx+1))+(draftarea1/2-totalareabelow)/length(jx);
  y1=(ylbottom1-y(idx-1));
  y2=(ylbottom1-y(idx));
  y3=(ylbottom1-y(idx+1));
  x1=x(idx-1);
  x2=x(idx);
  x3=x(idx+1);
  y(idx)=ylbottom1-(2*newarea-y1*(x2-x1)-y3*(x3-x2))./(x3-x1);
 end
end

plot(x,y,'b')

hold on

xxb=x;
yyb=y;

% calculate right half keel

xincrement=(scalefactorx*Lk2)./(2*maxx);
yincrement=(ylbottom2-ykbottom)./maxx;

x=zeros([1 maxx+1]);
y=zeros([1 maxx+1]);

while any(y(1:end-1)>=ylbottom2 | y(2:end)<=ykbottom)
 area=zeros([maxx 1]);

 for i=1:maxx+1
  x(i)=boxcenterxb+(i-1)*xincrement;
  coeff=(i/(maxx+1))^4;
  if i==1
   y(i)=ykbottom;
  elseif i==maxx+1
   y(i)=ylbottom2;
  else
   y(i)=y(i-1)+delt*yincrement*(rand-0.49);
   y(i)=max(ykbottom+3*yincrement,min(y(i),ylbottom2-3*yincrement));
   y(i)=(1-coeff)*y(i)+coeff*ylbottom2;
   area(i)=(x(i)-x(i-1))*(ylbottom2-0.5*(y(i)+y(i-1)));
  end
 end
 totalareabelow=sum(area(:));

 for idx=jx
  newarea=sum(area(idx:idx+1))+(draftarea2/2-totalareabelow)/length(jx);
  y1=(ylbottom2-y(idx-1));
  y2=(ylbottom2-y(idx));
  y3=(ylbottom2-y(idx+1));
  x1=x(idx-1);
  x2=x(idx);
  x3=x(idx+1);
  y(idx)=ylbottom2-(2*newarea-y1*(x2-x1)-y3*(x3-x2))./(x3-x1);
 end
end

plot(x,y,'b')

xxb=[xxb x(2:end)];
yyb=[yyb y(2:end)];


% calculate fill keel/ridge in with dots
ya=ykbottom:yincrementfirst:yrtop;
xk=zeros([length(xxb) length(ya)]);
yk=zeros([length(xxb) length(ya)]);
for i=1:length(xxb)
for j=1:length(ya)
 if ya(j)<=yyb(i)-yincrement | ya(j)>yrtop
  xk(i,j)=NaN;
  yk(i,j)=NaN;
 elseif ya(j)>=ylbottom1 & xxb(i)<(boxcenterxt-scalefactorx*Lr1/2)
  xk(i,j)=NaN;
  yk(i,j)=NaN;
 elseif ya(j)>=ylbottom2 & xxb(i)>(boxcenterxt+scalefactorx*Lr2/2)
  xk(i,j)=NaN;
  yk(i,j)=NaN;
 elseif mod(i-mod(j,2),2)==0
  xk(i,j)=NaN;
  yk(i,j)=NaN;
 else
  xk(i,j)=xxb(i);
  yk(i,j)=ya(j);
 end
end
end


% calculate left half ridge

xincrement=(boxcenterxt-(boxcenterxt-scalefactorx*Lr1/2))./maxx;
yincrement=(yrtop-ylsnow1)./maxx;

x=zeros([1 maxx+1]);
y=zeros([1 maxx+1]);

while any(y(2:end)<=ylsnow1 | y(1:end-1)>=yrtop)
 area=zeros([maxx 1]);
 for i=1:maxx+1
  x(i)=(boxcenterxt-scalefactorx*Lr1/2)+(i-1)*xincrement;
  coeff=(i/(maxx+1))^4;
  if i==1
   y(i)=ylsnow1;
  elseif i==maxx+1
   y(i)=yrtop;
  else
   y(i)=(y(i-1)+delt*yincrement*(rand-0.49));
   y(i)=max(ylsnow1+3*yincrement,min(y(i),yrtop-3*yincrement));
   y(i)=(1-coeff)*y(i)+coeff*yrtop;
   area(i)=(x(i)-x(i-1))*(0.5*(y(i)+y(i-1))-ylsnow1);
  end
 end
 totalareaabove=sum(area(:));

 for idx=jx
  newarea=sum(area(idx:idx+1))+((freeboardarea1/2)-totalareaabove)/length(jx);
  y1=(y(idx-1)-ylsnow1);
  y2=(y(idx)-ylsnow1);
  y3=(y(idx+1)-ylsnow1);
  x1=x(idx-1);
  x2=x(idx);
  x3=x(idx+1);
  y(idx)=ylsnow1+(2*newarea-y1*(x2-x1)-y3*(x3-x2))./(x3-x1);
 end
end

plot(x,y,'b')

xxt=[xlleft x];
yyt=[ylsnow1 y];


% calculate right half ridge

xincrement=(boxcenterxt-(boxcenterxt-scalefactorx*Lr2/2))./maxx;
yincrement=(yrtop-ylsnow1)./maxx;

x=zeros([1 maxx+1]);
y=zeros([1 maxx+1]);

while any(y(1:end-1)<=ylsnow2 | y(2:end)>=yrtop)
 area=zeros([maxx 1]);
 for i=1:maxx+1
  x(i)=boxcenterxt+(i-1)*xincrement;
  coeff=(i/(maxx+1))^4;
  if i==1
   y(i)=yrtop;
  elseif i==1
   y(i)=yrtop;
  else
   y(i)=(y(i-1)-delt*yincrement*(rand-0.49));
   y(i)=max(ylsnow2+3*yincrement,min(y(i),yrtop-3*yincrement));
   y(i)=(1-coeff)*y(i)+coeff*ylsnow2;
   area(i)=(x(i)-x(i-1))*(0.5*(y(i)+y(i-1))-ylsnow2);
  end
 end
 totalareaabove=sum(area(:));

 for idx=jx
  newarea=sum(area(idx:idx+1))+((freeboardarea2/2)-totalareaabove)/length(jx);
  y1=(y(idx-1)-ylsnow2);
  y2=(y(idx)-ylsnow2);
  y3=(y(idx+1)-ylsnow2);
  x1=x(idx-1);
  x2=x(idx);
  x3=x(idx+1);
  y(idx)=ylsnow2+(2*newarea-y1*(x2-x1)-y3*(x3-x2))./(x3-x1);
 end
end

plot(x,y,'b')

xxt=[xxt x xlright];
yyt=[yyt y ylsnow2];

x=[xxt fliplr(xxb)];
y=[yyt fliplr(yyb)];

% draw background patch of whitecol
minxx=min(x); maxxx=max(x);
minyy=min(y); maxyy=max(y);
patch([minxx minxx maxxx maxxx],[minyy maxyy maxyy minyy],whitecol,'EdgeColor',whitecol)
hold on

% draw shape of random ridge
patch(x,y,icecol,'EdgeColor',icecol)
hold on

% plot ridge line
x=[boxcenterxt boxcenterxt]; 
y=[yrtop ykbottom];
line(x,y,'Color',cols(2,:),'LineStyle',':');

% label ridge line
text(boxcenterxt,(yrtop+ykbottom)/2,'Ridge Line',...
           'Color',cols(2,:),...
           'VerticalAlignment','bottom',...
           'HorizontalAlignment','center',...
	   'Rotation',90,...
           'Interpreter','Latex',...
	   'BackgroundColor','none',...
           'Fontsize',labelsize-1,...
	   'Margin',0.5,...
           'EdgeColor','none')


% plot shear failure planes on right
seperation=0.025;

% downward slanting failure planes
thetar=atan2d((sealevely-ykbottom),-(xlright-boxcenterxb));
ystart=ykbottom+[-0.5:seperation:0.5];
xstart=xlright*ones(size(ystart));
for i=1:length(xstart)
 if ykbottom>ystart(i)
   radiusr(i)=max(0.,(ykbottom-ystart(i))/sind(thetar));
   xstart(i)=xstart(i)+radiusr(i).*cosd(thetar);
   ystart(i)=ykbottom;
 end
 if ystart(i)<yrtop & xstart(i)>boxcenterxb
  radiusr(i)=sqrt((yrtop-ystart(i)).^2+(boxcenterxb-xstart(i)).^2);
  xend(i)=xstart(i)+radiusr(i).*cosd(thetar);
  yend(i)=ystart(i)+radiusr(i).*sind(thetar);
  if yend(i)>yrtop
   radiusr(i)=(yrtop-ystart(i))/sind(thetar);
  elseif xend(i)<boxcenterxb
   radiusr(i)=(boxcenterxb-xstart(i))/cosd(thetar);
  end
  xend(i)=xstart(i)+radiusr(i).*cosd(thetar);
  yend(i)=ystart(i)+radiusr(i).*sind(thetar);
  line([xstart(i) xend(i)],[ystart(i) yend(i)],'Color',whitecol,'LineWidth',0.5)
 end
end

% upward slanting failure planes
thetar=atan2d((sealevely-ykbottom),(xlright-boxcenterxb));
ystart=ykbottom+[-0.5:seperation:0.5];
xstart=boxcenterxb*ones(size(ystart));
for i=1:length(xstart)
 if ykbottom>ystart(i)
   radiusr(i)=max(0.,(ykbottom-ystart(i))/sind(thetar));
   xstart(i)=xstart(i)+radiusr(i).*cosd(thetar);
   ystart(i)=ykbottom;
 end
 if ystart(i)<yrtop & xstart(i)<xlright
  radiusr(i)=sqrt((yrtop-ystart(i)).^2+(xlright-xstart(i)).^2);
  xend(i)=xstart(i)+radiusr(i).*cosd(thetar);
  yend(i)=ystart(i)+radiusr(i).*sind(thetar);
  if yend(i)>yrtop
   radiusr(i)=(yrtop-ystart(i))/sind(thetar);
  elseif xend(i)>xlright 
   radiusr(i)=(xlright-xstart(i))/cosd(thetar);
  end
  xend(i)=xstart(i)+radiusr(i).*cosd(thetar);
  yend(i)=ystart(i)+radiusr(i).*sind(thetar);
  line([xstart(i) xend(i)],[ystart(i) yend(i)],'Color',whitecol,'LineWidth',0.5)
  % label failure planes
  if ystart(i)>(ykbottom+seperation) & ystart(i)<=(ykbottom+2*seperation);
   text(0.99*(xstart(i)+xend(i))/2,0.99*(ystart(i)+yend(i))/2,'Failure Planes',...
        'FontSize',labelsize-1,'Color',whitecol,'HorizontalAlignment','center',...
        'VerticalAlignment','middle','BackgroundColor',icecol,'Margin',0.5,...
        'Rotation',thetar)
  end
 end
end

% label center line
centroid=false;
%centroid=true;
if centroid 
 ht=text(xc,sealevely-1.1*textoffset,'Centroid',...
	   'Interpreter','Latex',...
           'VerticalAlignment','top',...
           'HorizontalAlignment','center',...
           'Fontsize',labelsize-1,...
           'Color',cols(4,:));
 ext=get(ht,'Extent');
 mask=(xk>ext(1)-1.0*textoffset & xk<=ext(1)+ext(3)+1.0*textoffset & ...
       yk>=ext(2)-0.5*textoffset & yk<=(ext(2)+ext(4)-0.5*textoffset));
 xk(mask)=NaN;
 yk(mask)=NaN;
end

% plot sea level text across figure
st=text(xc-10*textoffset,sealevely-0.05*textoffset,...
           {'Sea Level'},...
	   'Interpreter','Latex',...
           'VerticalAlignment','top',...
           'HorizontalAlignment','right',...
           'Fontsize',labelsize-1,...
           'Color',bluecol,...
           'EdgeColor','none');
ext=get(st,'Extent');
mask=(xk>ext(1)-1.0*textoffset & xk<=ext(1)+ext(3)+1.0*textoffset & ...
      yk>=ext(2) & yk<=(ext(2)+ext(4)+1.0*textoffset));
xk(mask)=NaN;
yk(mask)=NaN;


% draw centroid
if centroid 
 plot(xc,yc,'o','Color',cols(4,:),'MarkerSize',10,'LineWidth',0.75)
 hp=plot(xc,yc,'+','Color',cols(4,:),'MarkerSize',20,'LineWidth',0.75)
end

drawnow


% plot feeding level ice plates left
hold on
x=[xlleft sealeft sealeft xlleft];
y=[ylsnow1 ylsnow1 ylbottom1 ylbottom1];
patch(x,y,icecol,'EdgeColor',icecol)
plot([sealeft sealeft],[ylsnow1 ylbottom1],':','Color',0.99999*[1 1 1],'LineWidth',3)

% plot feeding level ice plates right
hold on
x=[searight xlright xlright searight];
y=[ylsnow2 ylsnow2 ylbottom2 ylbottom2];
patch(x,y,icecol,'EdgeColor',icecol)
plot([searight searight],[ylsnow2 ylbottom2],':','Color',0.99999*[1 1 1],'LineWidth',3)

% plot feeding plate bottom line right
hold on
x=[xlright boxcenterxt+scalefactorx*Lr2/2];
y=[ylbottom2 ylbottom2];

% plot feeding plate bottom line left
hold on
x=[boxcenterxt-scalefactorx*Lr1/2 xlleft];
y=[ylbottom1 ylbottom1];

% plot feeding plate
%feeding=true;
feeding=false;

if feeding

 % add feeding plate left arrows
 x=[(boxcenterxt-scalefactorx*Lr1/2)-2*envelope ...
   (boxcenterxt-scalefactorx*Lr1/2)-2*arrowhead ...
   (boxcenterxt-scalefactorx*Lr1/2)-2*arrowhead ...
   (boxcenterxt-scalefactorx*Lr1/2) ...
   (boxcenterxt-scalefactorx*Lr1/2)-2*arrowhead ...
   (boxcenterxt-scalefactorx*Lr1/2)-2*arrowhead];
 y=[(ylbottom1+ylsnow1)/2 ...
   (ylbottom1+ylsnow1)/2 ...
   (ylbottom1+ylsnow1)/2+arrowhead ...
   (ylbottom1+ylsnow1)/2 ...
   (ylbottom1+ylsnow1)/2-arrowhead ...
   (ylbottom1+ylsnow1)/2];

 plot(x,y,'k')
 text((boxcenterxt-scalefactorx*Lr1/2)-2*envelope,y(1),...
      '$\dot{y}_{F_1}$',...
      'VerticalAlignment','middle',...
      'HorizontalAlignment','right',...
      'Interpreter','Latex',...
      'Fontsize',labelsize-1,...
      'Color','k')

 if rubblefloe
  text((xlleft+boxcenterxt-scalefactorx*Lr1/2)/2,y(1),...
      'floe ice',...
      'VerticalAlignment','middle',...
      'HorizontalAlignment','center',...
      'Interpreter','Latex',...
      'Fontsize',labelsize,...
      'Color','k')
 end

 % add feeding plate right arrows
 x=[(boxcenterxt+scalefactorx*Lr2/2)+2*envelope ...
   (boxcenterxt+scalefactorx*Lr2/2)+2*arrowhead ...
   (boxcenterxt+scalefactorx*Lr2/2)+2*arrowhead ...
   (boxcenterxt+scalefactorx*Lr2/2) ...
   (boxcenterxt+scalefactorx*Lr2/2)+2*arrowhead ...
   (boxcenterxt+scalefactorx*Lr2/2)+2*arrowhead];
 y=[(ylbottom2+ylsnow2)/2 ...
   (ylbottom2+ylsnow2)/2 ...
   (ylbottom2+ylsnow2)/2+arrowhead ...
   (ylbottom2+ylsnow2)/2 ...
   (ylbottom2+ylsnow2)/2-arrowhead ...
   (ylbottom2+ylsnow2)/2];

 plot(x,y,'k')
 text((boxcenterxt+scalefactorx*Lr2/2)+2*envelope+0.25*textoffset,y(1),...
      '$\dot{y}_{F_2}$',...
      'VerticalAlignment','middle',...
      'HorizontalAlignment','left',...
      'Interpreter','Latex',...
      'Fontsize',labelsize-1,...
      'Color','k')

 if rubblefloe
  text((xlright+boxcenterxt+scalefactorx*Lr2/2)/2,y(1),...
       'floe ice',...
       'VerticalAlignment','middle',...
       'HorizontalAlignment','center',...
       'Interpreter','Latex',...
       'Fontsize',labelsize,...
       'Color','k')
 end

end


% annotate angles of repose
radius=(scalefactorx*Lr1/2)*0.25;
theta=-[0:1:atan2d(ylbottom1-ykbottom,(boxrightx-boxleftx)/2)];
x=xlleft+radius*cosd(theta);
y=ylbottom1+radius*sind(theta);
line(x,y,'Color','k','LineWidth',0.5)

xi=[xlleft xlleft+7*textoffset];
yi=[ylbottom1 ylbottom1];
line(xi,yi,'Color','k','LineWidth',0.5)

alphatext1='$\alpha_{R}$';

text(xlleft+radius+0.5*textoffset,mean(y)-0.05*textoffset,...
           alphatext1,...
           'Color',greycol,...
           'VerticalAlignment','middle',...
	   'HorizontalAlignment','left',...
           'Interpreter','Latex',...
           'Fontsize',labelsize,...
           'EdgeColor','none')

theta=[0:1:atan2d(ylbottom1-ykbottom,(boxrightx-boxleftx)/2)];
x=boxcenterxt-scalefactorx*Lr1/2+radius*cosd(theta);
y=ylsnow1+radius*sind(theta);
line(x,y,'Color','k','LineWidth',0.5)

xi=[boxcenterxt-scalefactorx*Lr1/2 boxcenterxt-scalefactorx*Lr1/2+7*textoffset];
yi=[ylsnow1 ylsnow1];
line(xi,yi,'Color','k','LineWidth',0.5)

text(boxcenterxt-scalefactorx*Lr1/2+radius+0.5*textoffset,mean(y)+0.35*textoffset,...
           alphatext1,...
           'Color',greycol,...
           'VerticalAlignment','middle',...
           'HorizontalAlignment','left',...
           'Interpreter','Latex',...
           'Fontsize',labelsize,...
           'EdgeColor','none')

% plot sea level line
x=[sealeft searight]; 
y=[sealevely sealevely];
line(x,y,'Color',bluecol,'LineStyle','-.');

% label leading edge 
text(xlleft,(ylsnow1+ylbottom1)/2,'Leading Edge',...
           'Color',cols(2,:),...
           'VerticalAlignment','bottom',...
           'HorizontalAlignment','center',...
	   'Rotation',90,...
           'Interpreter','Latex',...
	   'BackgroundColor','none',...
           'Fontsize',labelsize-4,...
	   'Margin',1,...
           'EdgeColor','none')

% plot isostatic unit
x=[xlleft, xlleft, boxcenterxt-scalefactorx*Lr1/2, boxcenterxt,...
    boxcenterxt+scalefactorx*Lr2/2, xlright, xlright, boxcenterxb, xlleft];
y=[ylbottom1, ylsnow1, ylsnow1, yrtop, ylsnow2, ylsnow2, ylbottom2, ykbottom, ylbottom1];
plot(x,y,'-','Color',cols(2,:))
%plot([xlleft xlleft boxcenterxb xlright xlright boxcenterxt+scalefactorx*Lr2/2 boxcenterxt boxcenterxt-scalefactorx*Lr1/2],...
%     [ylsnow1 ylbottom1 ykbottom ylbottom2 ylsnow2 ylsnow2 yrtop ylsnow1],...
%     'o','Color',cols(2,:))
axis off
axis equal

hold on

% label isostatic unit
%isostaticlabel=true;
isostaticlabel=false;
if isostaticlabel
 text(xlleft+0.5*textoffset,ylsnow1,{'isostatic quanta'},...
	   'Interpreter','Latex',...
           'VerticalAlignment','bottom',...
           'HorizontalAlignment','left',...
           'Fontsize',labelsize,...
           'Color',cols(2,:))
end

if feeding
 % annotate ice feeding into ridge
 x=[xlleft-2*arrowhead xlleft+envelope ...
    xlleft+envelope xlleft+envelope+2*arrowhead ...
    xlleft+envelope xlleft+envelope];
 y=[(ylsnow1+ylbottom1)/2 (ylsnow1+ylbottom1)/2 ...
    (ylsnow1+ylbottom1)/2+arrowhead (ylsnow1+ylbottom1)/2 ...
    (ylsnow1+ylbottom1)/2-arrowhead (ylsnow1+ylbottom1)/2];
 plot(x,y,'Color',gcol)
 text(sealeft+0.05*textoffset,y(1),'$M_{F_1}$',...
     'Interpreter','Latex',...
     'VerticalAlignment','middle',...
     'Fontsize',labelsize-2,...
     'HorizontalAlignment','left',...
     'Color',gcol)

 x=[xlright+1.0*arrowhead xlright-envelope ...
    xlright-envelope xlright-envelope-2*arrowhead ...
    xlright-envelope xlright-envelope];
 y=[(ylsnow2+ylbottom2)/2 (ylsnow2+ylbottom2)/2 ...
    (ylsnow2+ylbottom2)/2+arrowhead (ylsnow2+ylbottom2)/2 ...
    (ylsnow2+ylbottom2)/2-arrowhead (ylsnow2+ylbottom2)/2];
 plot(x,y,'Color',gcol)
 text(searight-0.05*textoffset,y(1),'$M_{F_2}$',...
     'Interpreter','Latex',...
     'VerticalAlignment','middle',...
     'Fontsize',labelsize-2,...
     'HorizontalAlignment','right',...
     'Color',gcol)

 % annotate ice outflow in keel
 theta=[180:1:270]*pi/180;
 x=[boxcenterxb+envelope boxcenterxb+envelope ...
    boxcenterxb+envelope+2*envelope*(1+cos(theta)) ...
    boxcenterxb+envelope+2*envelope ...
    boxcenterxb+envelope+2*envelope boxcenterxb+envelope+2*envelope+2*arrowhead...
    boxcenterxb+envelope+2*envelope boxcenterxb+envelope+2*envelope];
 y=[ylbottom2-0.2*envelope (ylbottom2+ykbottom)/2 ...
    (ylbottom2+ykbottom)/2+2*envelope*(sin(theta)) ...
    (ylbottom2+ykbottom)/2-2*envelope ...
    (ylbottom2+ykbottom)/2-2*envelope-arrowhead (ylbottom2+ykbottom)/2-2*envelope ...
    (ylbottom2+ykbottom)/2-2*envelope+arrowhead (ylbottom2+ykbottom)/2-2*envelope];
 plot(x,y,'Color',gcol)
 
 theta=[0:-1:-90]*pi/180;
 x=[boxcenterxb-envelope boxcenterxb-envelope ...
    boxcenterxb-envelope-2*envelope*(1-cos(theta)) ...
    boxcenterxb-envelope-2*envelope ...
    boxcenterxb-envelope-2*envelope boxcenterxb-envelope-2*envelope-2*arrowhead...
    boxcenterxb-envelope-2*envelope boxcenterxb-envelope-2*envelope];
 y=[ylbottom1-0.2*envelope (ylbottom2+ykbottom)/2 ...
    (ylbottom2+ykbottom)/2+2*envelope*(sin(theta)) ...
    (ylbottom2+ykbottom)/2-2*envelope ...
    (ylbottom2+ykbottom)/2-2*envelope-arrowhead (ylbottom2+ykbottom)/2-2*envelope ...
    (ylbottom2+ykbottom)/2-2*envelope+arrowhead (ylbottom2+ykbottom)/2-2*envelope];
 plot(x,y,'Color',gcol)


 % annotate ice outflow in sail
 theta=[180:-1:90]*pi/180;
 x=[boxcenterxt+envelope boxcenterxb+envelope ...
    boxcenterxt+envelope+0.5*envelope*(1+cos(theta)) ...
    boxcenterxt+envelope+0.5*envelope ...
    boxcenterxt+envelope+0.5*envelope boxcenterxb+1.5*envelope+2*arrowhead...
    boxcenterxt+envelope+0.5*envelope boxcenterxb+1.5*envelope];
 y=[ylsnow2+0.2*envelope (2*ylsnow2+yrtop)/3 ...
    (2*ylsnow2+yrtop)/3+0.5*envelope*(sin(theta)) ...
    (2*ylsnow2+yrtop)/3+0.5*envelope ...
    (2*ylsnow2+yrtop)/3+0.5*envelope-arrowhead (2*ylsnow2+yrtop)/3+0.5*envelope ...
    (2*ylsnow2+yrtop)/3+0.5*envelope+arrowhead (2*ylsnow2+yrtop)/3+0.5*envelope];
 plot(x,y,'Color',gcol)

 theta=[0:1:90]*pi/180;
 x=[boxcenterxt-envelope boxcenterxb-envelope ...
    boxcenterxt-envelope-0.5*envelope*(1-cos(theta)) ...
    boxcenterxt-envelope-0.5*envelope ...   
    boxcenterxt-envelope-0.5*envelope boxcenterxb-1.5*envelope-2*arrowhead...
    boxcenterxt-envelope-0.5*envelope boxcenterxb-1.5*envelope];
 y=[ylsnow1+0.2*envelope (2*ylsnow2+yrtop)/3 ...
    (2*ylsnow2+yrtop)/3+0.5*envelope*(sin(theta)) ...
    (2*ylsnow2+yrtop)/3+0.5*envelope ...
    (2*ylsnow2+yrtop)/3+0.5*envelope-arrowhead (2*ylsnow2+yrtop)/3+0.5*envelope ...
    (2*ylsnow2+yrtop)/3+0.5*envelope+arrowhead (2*ylsnow2+yrtop)/3+0.5*envelope];
 plot(x,y,'Color',gcol)

end

% annotate lines next to angles

% plot line of frictional failure
x=[xlleft boxcenterxb];
y=[sealevely ykbottom];
plot(x,y,'Color','r','LineStyle','--')

% annotate friction angle
radius=(scalefactorx*Lr1/2)*0.5;
theta=-[0:1:atan2d(sealevely-ykbottom,(boxrightx-boxleftx)/2)];
x=xlleft+radius*cosd(theta);
y=sealevely+radius*sind(theta);
line(x,y,'Color','r','LineWidth',0.5)
text(xlleft+radius+0.5*textoffset,mean(y)-0.05*textoffset,...
           '$\psi_{R}$',...
           'Color','r',...
           'VerticalAlignment','middle',...
           'HorizontalAlignment','left',...
           'Interpreter','Latex',...
           'Fontsize',labelsize+2,...
           'EdgeColor','none')

% annotate horizontal friction
x=[boxcenterxb-2.2*envelope xlleft+7.5*envelope+2*arrowhead ...
   xlleft+7.5*envelope+2*arrowhead xlleft+7.5*envelope ...
   xlleft+7.5*envelope+2*arrowhead xlleft+7.5*envelope+2*arrowhead];
y=[(ylbottom1+ylsnow1)/2 (ylbottom1+ylsnow1)/2 ...
   (ylbottom1+ylsnow1)/2+arrowhead (ylbottom1+ylsnow1)/2 ...
   (ylbottom1+ylsnow1)/2-arrowhead (ylbottom1+ylsnow1)/2];
plot(x,y,'Color',cols(fricol,:))
text(mean([x(1) x(end)]),mean([y(1) y(end)])+0.4*textoffset,'$\sigma_{R_{\hat{x}}}$',...
     'Interpreter','Latex',...
     'VerticalAlignment','top',...
     'Fontsize',labelsize+4,...
     'Color',cols(fricol,:),...
     'HorizontalAlignment','center',...
     'EdgeColor','none')

% annotate vertical friction lower
x=[boxcenterxb-2.2*envelope ...
   boxcenterxb-2.2*envelope ...
   boxcenterxb-2.2*envelope+arrowhead ...
   boxcenterxb-2.2*envelope ...
   boxcenterxb-2.2*envelope-arrowhead ...
   boxcenterxb-2.2*envelope];

y=[ykbottom+2.2*envelope ...
   (ylbottom1+ylsnow1)/2-0.2*envelope-2*arrowhead ...
   (ylbottom1+ylsnow1)/2-0.2*envelope-2*arrowhead ...
   (ylbottom1+ylsnow1)/2-0.2*envelope ...
   (ylbottom1+ylsnow1)/2-0.2*envelope-2*arrowhead ...
   (ylbottom1+ylsnow1)/2-0.2*envelope-2*arrowhead];

plot(x,y,'Color',cols(fricol,:))

text(x(1),mean([y(1) y(2)]),'$\sigma_{R_{\hat{z}}}$',...
     'Interpreter','Latex',...
     'VerticalAlignment','bottom',...
     'Fontsize',labelsize+4,...
     'Rotation',90,...
     'Color',cols(fricol,:),...
     'HorizontalAlignment','center',...
     'EdgeColor','none')


info=false;
%info=true;
if info
 % plot information text in lower left corner
 info=true;
 if info 
  textbox={...
   ['$\phi_R=',num2str(porosity,'%8.1f'),'$'],...
   ['$\alpha_{R_1}=',num2str(alpha1),'^{\circ}, \alpha_{R_2}=',num2str(alpha2),'^{\circ}$'],...
   ['$h_{F_{i_1}}=',num2str(hfi1,'%8.1f'),'$ m, $h_{F_{i_2}}=',num2str(hfi2,'%8.1f'),'$ m'],...
   ['$h_{F_{s_1}}=',num2str(hfs1,'%8.1f'),'$ m, $h_{F_{s_2}}=',num2str(hfs2,'%8.1f'),'$ m'],...
   ['$H_K=',num2str(Hk,'%8.1f'),'$ m, $H_S=',num2str(Hr,'%8.1f'),'$ m']
   };
 end

 if info
  text(sealeft+textoffset,ykbottom+0.0*textoffset,textbox,...
	    'Interpreter','Latex',...
            'VerticalAlignment','bottom',...
            'HorizontalAlignment','left',...
            'Fontsize',labelsize-2,...
            'Color','k',...
            'EdgeColor','none')
 end

 % provide horizontal and vertical scale bar in lower right corner
 x=[searight-2*scalefactorx searight]-5*textoffset;
 y=[ykbottom ykbottom]+5*textoffset;
 line(x,y,'Color',0.0*[1 1 1],'LineStyle','-');
 line([x(1) x(1)],[y(1) y(1)+arrowhead],...
           'Color',0.0*[1 1 1],'LineStyle','-');
 text(sum(x)/2,sum(y)/2-0.5*textoffset,...
           ['2 m'],...
	   'Interpreter','Latex',...
           'VerticalAlignment','bottom',...
           'Fontsize',labelsize-1,...
           'HorizontalAlignment','center',...
           'EdgeColor','none')

 x=[searight searight]-5*textoffset;
 y=[ykbottom ykbottom+2*scalefactory]+5*textoffset;
 line(x,y,'Color',0.0*[1 1 1],'LineStyle','-');
 line([x(1)-arrowhead x(1)],[y(2) y(2)],...
          'Color',0.0*[1 1 1],'LineStyle','-');
 text(sum(x)/2+0.5*textoffset,sum(y)/2,...
           ['2 m'],...
	   'Interpreter','Latex',...
           'VerticalAlignment','bottom',...
           'HorizontalAlignment','center',...
           'Fontsize',labelsize-1,...
           'Rotation',90,...
           'EdgeColor','none')

end

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



