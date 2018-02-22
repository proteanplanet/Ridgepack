function [AX,H2]=ridgepack_plotpsd(faxis,psd,leg,conf,cutconf,description,units,flo,fhi,plo,phi,col,pattern)

% ridgepack_plotpsd - Plot one- or two-sided power spectral density
%
% function [AX,H2]=ridgepack_plotpsd(faxis,psd,leg,conf,cutconf,description,units,flo,fhi,plo,phi,col,pattern)
% 
% This function plots one-sided or two-sided (rotary) spectra given the following inputs:
%
% Input:
% faxis       - Frequency axis with size (Nx1). If faxis has values less than zero,
%               a rotary spectra is plotted, otherwise a scalar spectra is plotted.
% psd         - Power spectral density corresponding to each faxis value (NxM)
%               for a total of M traces.
% leg         - Cell array of text for a legend (enter an empty vector [] for no legend)
%               Only include legend entries for the lines you wish to include in the 
%               legend, and for the others, include a null character (e.g.{'t1','','t2')
%               which will only plot the first and third lines in the legend.
% conf        - Lower and upper confidence bounds on psd with size (Nx2) 
%               (the first index of the second dimension is the lower error bound, 
%               the second index of second dimension is the upper error bound)
% cutconf     - Confidence limit for confidence bounds in conf
% description - A cell array of strings describing spectra details (optional)
% units       - Units of psd.  If this is omitted, psd is plotted in dB, otherwise
%               it is plotted on semilogarithmic axes. (optional)
% flo         - Minimum frequency plotted (optional)
% fhi         - Maximum frequency plotted (optional)
% plo         - Minimum psd axis value in dB/Bandwidth (optional)
% phi         - Maximum psd axis value in dB/Bandwidth (optional)
% col         - Color of the trace(s) to be plotted in cell array (e.g.{[0 0 1],[0 1 ]})
% pattern     - Pattern required for trace in cell array (e.g. {'--',':'})
%
% Output:
% AX - axis handles
% H2 - psd handles
%
% Andrew Roberts, Naval Postgraduate School, March 2018  (afrobert@nps.edu)
%

global debug;
if debug; disp(['Entering ',mfilename,'...']); end

if ~ishold; clf; end

% check necessary inputs
if nargin<2
 error('Must supply faxis and psd at the very least')
end

if nargin==4
 error('Must specify the confidence level with the conf array')
end


% determine if this is rotary and wrap faxis if it is
if min(faxis)<0; 
 disp('Plotting rotary (two sided) spectra'); 
 faxis=cat(1,-faxis(end),faxis);
 psd=cat(1,psd(end,:),psd(:,:));
 if nargin>3; conf=cat(1,conf(end,:),conf); end
end

% set axis limits
if nargin<9
 %flo=ceil(min(faxis));
 %fhi=floor(max(faxis));
 flo=max(-4.5,ceil(min(faxis)));
 fhi=min(4.5,floor(max(faxis)));
end
if nargin<11
 plo=5*floor(min(ncdb(psd(faxis>flo & faxis<fhi)))/5);
 phi=10*ceil(max(ncdb(psd(faxis>flo & faxis<fhi)))/10);
end

% plot graph initially on a linear scale with units of dB
psd=ncdb(psd); 

% plot one side in dB and other in semilogy. This is done this way
% so that text and the error bar can be plotted more easily on the 
% linear scale (with units of dB).
[AX,H1,H2]=plotyy(faxis,psd,faxis,ncinvdb(psd),'plot','semilogy'); hold on;
if nargin>11

 for hn=1:length(H1)
  set(H1(hn),'Color',col{hn}); set(H2(hn),'Color',col{hn}); 
 end

 if nargin>12
  for hn=1:length(H1)
   set(H1(hn),'LineStyle',pattern{hn}); set(H2(hn),'LineStyle',pattern{hn}); 
  end
 end

 % create single group for the legend (only using second handle H2)
 hSGroup=hggroup; 
 set(H2,'Parent',hSGroup)
 set(hSGroup,'Visible','off')

 hCGroup=hggroup; 
 set(H1,'Parent',hCGroup)
 set(get(get(hCGroup,'Annotation'),'LegendInformation'),'IconDisplayStyle','on');

elseif size(psd,2)==1 
 set(H1,'Color','b'); set(H2,'Color','b'); 
end

% create legend only for lines with non-empty legend entries
if  ~isempty(leg)
 legcount=0;
 for pleg=1:length(leg)
  if ~isempty(leg{pleg})
   legcount=legcount+1;
   hleg(legcount)=H2(pleg);
   tleg{legcount}=leg{pleg};
  end
 end
 legend_handle=legend(hleg,tleg);
 set(legend_handle, 'Box', 'off');
 set(legend_handle, 'Color', 'none');
end

xlim(AX(1),[flo, fhi]); xlim(AX(2),[flo, fhi]);
ylim(AX(1),[plo, phi]); ylim(AX(2),[ncinvdb(plo), ncinvdb(phi)]);
xlab=xlabel(['Frequency (cycles day$^{-1}$)']);

% If units are not defined, use a relative axis with dB 
if nargin<7 | isempty(units) | strcmp(units,'')
 disp('Using linear dB axis')
 set(AX(1),'YAxisLocation','left','YColor','k','Ytick',[-100:5:100]);
 set(AX(2),'YAxisLocation','right','YColor','k','YTickLabel',[],'YMinorTick','off');
 ylabel(AX(1),'Power Spectral Density (dB)')
else % switch the axes to semilogy
 disp('Using semilog axis')
 set(AX(2),'YAxisLocation','left','YColor','k','Ytick',10.^[-20:1:3]);
 set(AX(1),'YAxisLocation','right','YColor','k','YTickLabel',[],'Ytick',[]);
 ylabel(AX(2),['Power Spectral Density (',units,')'])
end

position=get(AX(1),'Position');
fontsize=min(max(8,12*position(4).^(1/3)),10);


% if rotary, add clockwise and anticlockwise indicators
if min(faxis)<0
 pos=get(xlab,'position');
 T(1)=text(flo,pos(2),'$\hookleftarrow$ Clockwise',...
      'HorizontalAlignment','left','VerticalAlignment','cap','Parent',AX(1));
 T(2)=text(fhi,pos(2),'Anticlockwise $\hookrightarrow$',...
      'HorizontalAlignment','right','VerticalAlignment','cap','Parent',AX(1));
end

% provide error bar
if nargin>4
 disp('Plotting error bar')
 conf=ncdb(conf*ncinvdb(psd(1,1)));
 errbelow=psd(1,1)-conf(1,1);
 errabove=conf(1,2)-psd(1,1);
 if flo<0; xerror=flo+(fhi-flo)*0.055; else; xerror=flo+(fhi-flo)*0.50; end
 %yerror=plo+(phi-plo)*0.78; 
 yerror=plo+(phi-plo)*0.65; 
 hCLines=errorbar(xerror,yerror,errbelow,errabove,'.k'); hold on;
 text(xerror,yerror,['$\;$ ',num2str(cutconf),'\% confidence'])

 % remove error bar from legend information
 hCGroup=hggroup; 
 set(hCLines,'Parent',hCGroup)
 set(get(get(hCGroup,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
end;

% provide processing information
if nargin>=6 
 disp('Adding description')
 N=size(psd,1);
 if flo<0; xinfo=flo+(fhi-flo)*0.05; else; xinfo=flo+(fhi-flo)*0.10; end
 yinfo=plo+(phi-plo)*0.90;
 text(xinfo,yinfo,description,'HorizontalAlignment','left');
end

% remove hold on plotting
hold off;

set(AX(1),'FontSize',fontsize,'box','on');
set(AX(2),'FontSize',fontsize,'box','on');

drawnow

if debug; disp(['...Leaving ',mfilename]); end


