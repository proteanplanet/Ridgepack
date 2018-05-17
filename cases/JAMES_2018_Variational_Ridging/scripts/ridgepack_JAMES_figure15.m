% ridgepack_JAMES_figure15 - Generates Figure 15 in JAMES Variation Ridging paper
% 
% This script generates Figure 15 from:
%
% Roberts, A.F., E.C. Hunke, S.M. Kamal, W.H. Lipscomb, C. Horvat, W. Maslowski (2018),
% Variational Method for Sea Ice Ridging in Earth System Models, Part I: Theory, 
% submitted to J. Adv. Model Earth Sy.
%
% This function generates and plot a zeta-hat plane as a test of the function 
% ridgepack_zetahatplane in the morphology library.  The script also calculates
% mean values of porosity, strain, and angle of repose for the entire model
% and with observational limitations imposed by Peter Wadhams when interpreting
% sea ice draft measurements from submarine SONAR.
%
% Andrew Roberts, Naval Postgraduate School, March 2018 (afrobert@nps.edu)

clear
close all

% set ice and snow thickness
hfi=[0.5 2.0];
hfs=[0 0];

% color along line
cont=[0:0.05:0.40];
cmap=ridgepack_colormap(cont,0,'parula');
cmap=flipud(cmap);
colormap(cmap)

% Limits and labels
xlab='Ridge Width, $L_K$ (m)';
xmin=10^0;
xmax=2*10^3;
ylab='Probability Density';
ymin=10^-9;
ymax=10^-2;

% switch to only indicate gradient once 
final=true;
%final=false;

for i=1:length(hfi)

 % calculate path for a given thickness
 [strainp,phip,alphap,VR,HK,HS,LK,LS]=ridgepack_trajectory(hfi(i),0);

 % work per ridge shape
 energyratio=sum(LK(1:end).*VR(1:end))./(LK(1:end).*VR(1:end));

 % probability
 probability=energyratio./sum(energyratio);

 % explanatory output
 disp(['------- For floe ice thickness ',num2str(hfi(i),'%8.1f'),...
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
 ridgepack_text(x(y<2*10^-8),y(y<2*10^-8),['$h_F{=}',num2str(hfi(i),...
                 '%5.1f'),'$m'],9,cmap(end,:));

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

ridgepack_colorbar(cont,'\phi_R');

% determine directory for read/write
dir=fileparts(which(mfilename));
cd([dir(1:strfind(dir,'scripts')-1),'output']);

% determine filename
x=strfind(mfilename,'_');
thisfilename=mfilename;
graphicsout=thisfilename(x(end)+1:end);

% output
disp(['Writing graphics output ',graphicsout,' to:',char(13),' ',pwd])

% print figure
ridgepack_fprint('epsc',graphicsout,1,2)


