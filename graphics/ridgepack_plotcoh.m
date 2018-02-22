function ridgepack_plotcoh(faxis,coh,cohou,leg,conf,cutconf,description,flo,fhi)

% ridgepack_plotCOH - Plots coherence squared spectra
%
% function ridgepack_plotcoh(faxis,coh,cohou,leg,conf,cutconf,description,flo,fhi)
%
% This function plots coherence spectra for standard spectral 
% coherence squared, rotary coherence or auto coherence given 
% information to plot M traces:
%
% Input:
% faxis       - frequency axis with size (Nx1)
% coh         - co-rotating coherence squared spectra (NxM)
% cohou       - counter-rotating coherence squared spectra (NxM)
%               This can be entered as an empty array [] if 
%               rotary coherence is not being plotted or if you only
%               wish to plot inner coherence of the rotary component.
% leg         - cell array of legend text for each trace (empty for one trace)
% conf        - confidence array (NxM) (optional)
% cutconf     - confidence array value (e.g. 95% or 99%) (optional)
% description - a cell array of strings describing spectra details (optional)
% flo         - minimum frequency plotted (optional)
% fhi         - maximum frequency plotted (optional)
%
% Output is graphical only.
%
% Andrew Roberts, Naval Postgraduate School, March 2018  (afrobert@nps.edu)

global debug;
if debug; disp(['Entering ',mfilename,'...']); end

% check inputs
if nargin<2
 error('Must supply faxis and coh at the very least')
elseif nargin<3
 cohou=[];
end

if nargin<5
 conf=zeros(size(coh));
elseif nargin==5
 error('Must specify the confidence level with the conf array')
end

% determine if this is rotary and wrap faxis if it is
if min(faxis)<0;
 disp('Plotting rotary (two sided) spectra');
 faxis=cat(1,-faxis(end),faxis);
 coh=cat(1,coh(end,:),coh(:,:));
 if ~isempty(cohou); cohou=cat(1,cohou(end,:),cohou(:,:)); end
 conf=cat(1,conf(end),conf);
end

if nargin<7
 description={' '};
end

%if nargin<8; flo=ceil(min(faxis)); end
%if nargin<9; fhi=floor(max(faxis)); end
if nargin<8; flo=-4.5; end
if nargin<9; fhi=4.5; end

if isempty(cohou)

 disp('Not plotting outer coherence');

 plot(faxis,conf,':k') ; hold on;
 plot(faxis,coh,'LineWidth',0.11); hold off; 
 if flo<0 & fhi>0 
  ylabel('Inner coherence')
 else
  ylabel('Coherence')
 end
 xlab=xlabel('Frequency (cycles {day$^{-1}$})');
 axis([flo fhi 0 1.05])
 
 text(flo+(fhi-flo)*0.02,1.03,description,'VerticalAlignment','top')

 if flo<0 & fhi>0 
  pos=get(xlab,'position');
  text(flo,pos(2),'$\leftarrow$ clockwise',...
      'HorizontalAlignment','left','VerticalAlignment','cap')
  text(fhi,pos(2),'anticlockwise $\rightarrow$',...
      'HorizontalAlignment','right','VerticalAlignment','cap')
 end

else

 a1=subplot('Position',[0.10,0.53,0.80,0.40]); plot(faxis,conf,':k') ; hold on;
 plot(faxis,coh); hold off; 
 ylabel('Inner coherence')
 set(gca,'XtickLabel',[]);

 flo 
 fhi
 axis([flo fhi 0 1.0])

 a2=subplot('Position',[0.10,0.10,0.80,0.40]); plot(faxis,conf,':k') ; hold on;
 plot(faxis,cohou); hold off;
 xlab=xlabel('Frequency (cycles {day$^{-1}$})');
 ylabel('Outer coherence')

 linkaxes([a1 a2], 'x')

 axis([flo fhi 0 1.0])

 text(flo+(fhi-flo)*0.04,0.90,description,'VerticalAlignment','top')

 pos=get(xlab,'position');
 text(flo,pos(2),{' ',' ','$\leftarrow$ clockwise'},...
      'HorizontalAlignment','left','VerticalAlignment','top','Parent',a2)
 text(fhi,pos(2),{' ',' ','anticlockwise $\rightarrow$'},...
      'HorizontalAlignment','right','VerticalAlignment','top','Parent',a2)

end

if nargin>5 & size(coh,2)>1 ; 
        leg={[num2str(cutconf),'% confidence'],leg{1:end}};
	legend(leg) ; legend('boxoff') ; 
elseif nargin>5 ;
	legend([num2str(cutconf),'% confidence']) ; legend('boxoff') ; 
elseif nargin>3 & size(coh,2)>1 ;
	legend(leg) ; legend('boxoff') ; 
end

% set the current axes to the top graph if plotting inner and outer coherence
if ~isempty(cohou); set(gcf,'CurrentAxes',a1); end

drawnow

if debug; disp(['...Leaving ',mfilename]); end


