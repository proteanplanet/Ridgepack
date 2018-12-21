% ridgepack_JAMES_figure15 - Generates Figure 15 in JAMES Variational Ridging paper
% 
% This script generates Figure 15 and Table 2 from:
%
% Roberts, A.F., E.C. Hunke, S.M. Kamal, W.H. Lipscomb, C. Horvat, 
% W. Maslowski (2019), Variational Method for Sea Ice Ridging in 
% Earth System Models, J. Adv. Model Earth Sy.
%
% This function generates and plot a zeta-hat plane as a test of the function 
% ridgepack_zetahatplane in the morphology library.  The script also calculates
% mean values of porosity, strain, and angle of repose for the entire model
% and with observational limitations imposed by Peter Wadhams when interpreting
% sea ice draft measurements from submarine SONAR.
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

% porosity color divisions along line
cont=[0:0.05:0.40];
cmap=ridgepack_colormap(cont,0,'parula');
cmap=flipud(cmap);
colormap(cmap)

% set ice and snow thickness
hf=[0.5 2.0];
hfs=[0 0];

% Limits and labels
xlab='Ridge Width, $\hat{L}_K$ (m)';
xmin=10^0;
xmax=7*10^3;
ylab='Probability Density';
ymin=2*10^-11;
ymax=10^-2;

% switch on to only indicate gradient once, rather than for both lines
final=true;
%final=false;

% plot probabilities of ridges forming
for i=1:length(hf)

 % calculate path for a given thickness
 [strainp,phip,alphap,VR,HK,HS,LK,LS]=ridgepack_trajectory(hf(i),0,0.001);

 % work per ridge shape
 energyratio=sum(LK(1:end).*VR(1:end))./(LK(1:end).*VR(1:end));

 % probability
 probability=energyratio./sum(energyratio);

 % explanatory output
 disp(['------- For floe ice thickness ',num2str(hf(i),'%8.1f'),...
       'm ----------------------'])

 % calculate mean porosity
 disp(['Mean porosity: ',num2str(sum(probability.*phip),'%8.2f')])

 % calculate mean strain
 disp(['Mean strain: ',num2str(sum(probability.*strainp),'%8.2f')])

 % calculate mean alpha
 disp(['Mean angle of repose: ',num2str(sum(probability.*alphap),'%8.2f')])

 % calculate 5m cutoff
 mask=(HK>5.0);
 maskprob=(probability.*mask)/sum(probability.*mask);

 % calculate mean porosity
 disp(['Mean porosity with Wadhams 5m draft cutoff: ',num2str(sum(maskprob.*phip),'%8.2f')])

 % calculate mean strain
 disp(['Mean strain with Wadhams 5m draft cutoff: ',num2str(sum(maskprob.*strainp),'%8.2f')])

 % calculate mean alpha
 disp(['Mean angle of repose with Wadhams 5m draft cutoff: ',num2str(sum(maskprob.*alphap),'%8.2f')])

 % x coordinate is total ridge thickness
 x=LK;

 % y coordinates is PDF
 y=probability;

 % cumulative distribution
 cumulativedist=cumsum(probability,'forward');

 %plot(x,y,'bo')
 x2=log(x);
 y2=log(y);
 p2=polyfit(x2,y2,1);
 f2=polyval(p2,x2);
 d2=(5-p2(1))/2;
 sstot=sum((y2-mean(y2)).^2);
 ssres=sum((y2-f2).^2);
 rsquared2 = 1-(ssres/sstot);
 hh(i)=plot(exp(x2),exp(f2),':r');
 legtext{i}=['D=',num2str(p2(1),'%5.2f'),', R^2=',num2str(rsquared2,'%5.3f')];

 if i==1
  set(gca,'Xscale','log')
  set(gca,'Yscale','log')
  xlabel(xlab)
  ylabel(ylab)
  xlim([xmin xmax])
  ylim([ymin ymax])
  hold on
 end

 z = zeros(size(x));
 col = phip;  % This is the color, vary with x in this case.
 surface([x;x],[y;y],[z;z],[col;col],...
        'facecol','no',...
        'edgecol','interp',...
        'linew',2);
 ridgepack_text(x(y<2*10^-9),y(y<2*10^-9),['$h_F{=}',num2str(hf(i),...
                 '%5.1f'),'$m'],10,cmap(end,:));

end

% tidy up written output
disp(['----------------------------------------------------------'])

% only plot first line data because they both agree
if final
 legend(hh(1),legtext{1},'Location','northeast')
 legend('boxoff')
else
 legend(hh,legtext,'Location','northeast')
 legend('boxoff')
end

ridgepack_colorbar(cont,'\phi_R')

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

% Calculate mean values for table
hf=[0.2 0.5 1.0 2.0 5.0 5.5];

for i=1:length(hf)

 % calculate path for a given thickness
 [EPSILON,PHI,ALPHAHAT,VR,HK,HS,LK,LS]=...
     ridgepack_trajectory(hf(i),0,0.001,0.01,0.99);

 % work per ridge shape
 energyratio=sum(LK(:).*VR(:))./(LK(:).*VR(:));

 % probability
 probability=energyratio./sum(energyratio);

 % calculate mean porosity
 meanphi=sum(probability.*PHI(:));

 % calculate mean alpha
 meanalpha=sum(probability.*ALPHAHAT(:));

 % calculate 5m cutoff
 mask=(HK(:)>5.0);
 maskprob=(probability.*mask)/sum(probability.*mask);

 % calculate mean porosity
 meanphilim(i)=sum(maskprob.*PHI(:));

 % calculate mean alpha
 meanalphalim(i)=sum(maskprob.*ALPHAHAT(:));

end

disp(['----------------------------------------------------------'])
disp('Information for published Table 2 of Roberts et al. (2019)')
disp(['$h_F$ (m) ',num2str(hf,'& %6.2f '),' \\'])
disp(['$\alpha_{\mathrm{obs}}$ ',num2str(meanalphalim,'& %6.1f '),' \\'])
disp(['$\phi_{\mathrm{obs}}$  ',num2str(meanphilim,'& %6.2f '),' \\'])
disp(['----------------------------------------------------------'])

