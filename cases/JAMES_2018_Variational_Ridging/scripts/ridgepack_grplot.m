% ridgepack_grplot - Generates plot of gR on discrete and continuous grid
% 
% This function plots a comparison of the step function building block of 
% redistribution on both a discrete and continuous grid as a support for 
% checking numerics used in the paper:
%
% Roberts, A.F., E.C. Hunke, S.M. Kamal, W.H. Lipscomb, C. Horvat, W. Maslowski (2018),
% Variational Method for Sea Ice Ridging in Earth System Models, Part I: Theory, 
% submitted to J. Adv. Model Earth Sy.
%
% Andrew Roberts, Naval Postgraduate School, April 2018 (afrobert@nps.edu)

clear
close all

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

% determine directory for read/write
dir=fileparts(which(mfilename));
outdir=[dir(1:strfind(dir,'scripts')-1),'output'];
[status,msg]=mkdir(outdir);
cd(outdir);

% determine filename
graphicsout=mfilename;

% output
disp(['Writing graphics output ',graphicsout,' to:',char(13),' ',pwd])

% print figure
ridgepack_fprint('epsc',graphicsout,1,2)
ridgepack_fprint('png',graphicsout,1,2)

