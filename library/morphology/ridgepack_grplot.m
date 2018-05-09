function ridgepack_grplot

% RIDGEPACK_GRPLOT - Plot gR on discrete and continuous grid
% 
% function ridgepack_grplot
% 
% This function plots a comparison of the step function building 
% block of redistribution on a discrete and continuous grid.
%
% Ridgepack Version 1.0.
% Andrew Roberts, Naval Postgraduate School, March 2018 (afrobert@nps.edu)

hF=2;
hFs=0.3;
epsilon=-1/10;
phi=0.2;

% Calculate discrete thickness distribution on given thickness axis
[hincr,eincr,hgrid]=ridgepack_gridinit;
[GRHPHI]=ridgepack_grhphi(hF,hFs,epsilon,phi);

% check area under the step function
disp(['Total area of gR on grid is: ',num2str(sum(hincr.*GRHPHI))]);

% plot function
clf
plot(hgrid,GRHPHI)
hold on

% Thickness distribution in continuous space
[VR,ALPHAHAT,HK,HS,LK,LS]=ridgepack_energetics(hF,hFs,epsilon,phi);

step(1)=hF;
gR(1)=0;

step(2)=hF;
gR(2)=2/(LK*tand(ALPHAHAT));

step(3)=hF+(LK-LS)*tand(ALPHAHAT)/2;
gR(3)=2/(LK*tand(ALPHAHAT));

step(4)=hF+(LK-LS)*tand(ALPHAHAT)/2;
gR(4)=1/(LK*tand(ALPHAHAT));

step(5)=hF+(LK+LS)*tand(ALPHAHAT)/2;
gR(5)=1/(LK*tand(ALPHAHAT));

step(6)=hF+(LK+LS)*tand(ALPHAHAT)/2;
gR(6)=0;

% check area under the step function
area=gR(3)*(LK-LS)*tand(ALPHAHAT)/2+gR(5)*(2*LS)*tand(ALPHAHAT)/2;
disp(['Total area of continuous gR is: ',num2str(area)]);

% plot function
plot(step,gR)
legend({'Discrete','Continuous'})
xlabel('$h$ (m)')
ylabel('$g_{R}(h,\phi)$')
title(['$\epsilon_{R_I}$=',num2str(epsilon),', $\phi$=',num2str(phi)]);
xlim([max(0,step(1)-0.5) step(6)+0.5])

% determine directory for read/write of zeta-hat plane data
dir=fileparts(which(mfilename));
cd([dir(1:strfind(dir,'library')-1),'cases/JAMES_2018_Variational_Ridging/output'])
disp(['Writing graphics output to ',pwd])

x=strfind(mfilename,'_')
p=x(end)+1
xxx=mfilename;
xxx(p:end)
%graphicsout=mfilename(x(end)+1:end)


% print figure
ridgepack_fprint('png','ridgepack_grplot',1,2)
ridgepack_fprint('epsc','ridgepack_grplot',1,2)

