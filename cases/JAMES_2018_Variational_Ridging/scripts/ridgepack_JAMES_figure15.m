
clear
close all

% set ice and snow thickness
hfi=[0.5 2.0];
hfs=0;

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
 [hfii,strainp,pep,phip,alphap,HK,HS,LK,LS]=paper_ridge_expected_plane(true,hfi(i));

 % work per ridge shape
 energyratio=sum(LK(1:end).*pep(1:end))./(LK(1:end).*pep(1:end));

 % probability
 probability=energyratio./sum(energyratio);

 % calculate mean porosity
 sum(probability.*phip)

 % calculate mean strain
 sum(probability.*strainp')

 % calculate mean alpha
 sum(probability.*alphap)

 % calculate 5m cutoff
 mask=(HK>5.0);
 maskprob=(probability.*mask)/sum(probability.*mask);

 % calculate mean porosity
 sum(maskprob.*phip)

 % calculate mean strain
 sum(maskprob.*strainp')

 % calculate mean alpha
 sum(maskprob.*alphap)

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
 hh(i)=plot(exp(x2),exp(f2),'--b')
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
 ridgepack_text(x(y<2*10^-8),y(y<2*10^-8),['$h_F{=}',num2str(hfi(i),'%5.1f'),'$m'],9,cmap(end,:))

end

% only plot first line data because they both agree
if final
 legend(hh(1),legtext{1},'Location','northeast')
 legend('boxoff')
else
 legend(hh,legtext,'Location','northeast')
 legend('boxoff')
end


ridgepack_colorbar(cont,'\phi_R')

cd /Users/aroberts/Publications/2015_Unified_Morphology_1/figures 

ridgepack_fprint('png',['Ridge_Formation_Probability_lowres.png'],1,1)
ridgepack_fprint('epsc',['Ridge_Formation_Probability_lowres.eps'],1,1)
ridgepack_fprint('png',['Ridge_Formation_Probability.png'],1,2)
ridgepack_fprint('epsc',['Ridge_Formation_Probability.eps'],1,2)

