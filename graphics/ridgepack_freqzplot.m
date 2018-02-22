function ridgepack_freqzplot(H,W,b,a,description,units)

% ridgepack_freqzplot - Plots the magnitude and phase of a given frequency response
%
% function ridgepack_freqzplot(H,W,b,a,description,units)
% 
% Input:
% H           - frequency response 
% W           - frequency in cycles per units
% b           - b coefficient in difference equation
% a           - a coefficient in difference equation
% description - information for the title of the plot. 
% units       - time or space units of series (days, months, km etc)
%
% Output:
% Output is graphical.
%
% Andrew Roberts, Naval Postgraduate School, March 2018  (afrobert@nps.edu)
%

global debug;
if debug; disp(['Entering ',mfilename,'...']); end

clf

maxy=max(1.05,max(abs(H))); % maximum abs(H);

subplot('Position',[0.10,0.70,0.80,0.17]), plot(W,abs(H));
ylim([0 maxy]);
xlim([min(W),max(W)]);
set(gca,'XtickLabel',[]);
ylabel('Magnitude (linear)');
grid on;

if isinf(maxy) | isnan(maxy)
 title('WARNING: Magnitude has Infinite or NaN values')
else
 title(description,'HorizontalAlignment','left','position',[0.0 maxy*1.1/1.05 1])
end

subplot('Position',[0.10,0.50,0.80,0.17]), [AX,H1,H2]=plotyy(W,ncdb(abs(H)),W,abs(H),'plot','semilogy');
ylim(AX(1),[ncdb(1.e-5),ncdb(maxy)]);
ylim(AX(2),[1.e-5,maxy]);
xlim(AX(1),[min(W),max(W)]);
xlim(AX(2),[min(W),max(W)]);
set(AX(1),'XtickLabel',[],'Ytick',[-130:10:0],'YColor','k');
set(AX(2),'XtickLabel',[],'Ytick',[1.e-10 1.e-9 1.e-8 1.e-7 1.e-6 1.e-5 1.e-4 1.e-3 1.e-2 1.e-1 1]);
set(AX(2),'YColor','k','YMinorTick','off');
set(H1,'Color','b');
set(H2,'Color','b');
ylabel(AX(1),'|H({\theta})| (dB)');
ylabel(AX(2),'|H({\theta})| (log)');
grid on;

subplot('Position',[0.10,0.35,0.80,0.12]), plot(W,angle(H)*180/pi);
xlim([min(W),max(W)]);
xlabel(['Frequency (',units,'^{-1})']);
ylabel('{\angle}H({\theta}) (degrees)');
grid on;

if a==1; 
	subplot('Position',[0.10,0.07,0.80,0.15]), stem(b,'b.');
	ylabel('h(n)'); 
	xlabel('n');
	xlim([0,length(b)]);
	ylim([min(b)-max(b)/10,max(b)+max(b)/10]);
else; 
	subplot('Position',[0.10,0.07,0.35,0.15]), stem(b,'b.');
	ylabel('b(n)'); 
	xlabel('n');
	xlim([0,length(b)]);
	ylim([min(b)-max(b)/10,max(b)+max(b)/10]);

	subplot('Position',[0.55,0.07,0.35,0.15]), stem(a,'b.');
	ylabel('a(n)'); 
	xlabel('n');
	xlim([0,length(a)]);
	ylim([min(a)-max(a)/10,max(a)+max(a)/10]);
end
grid on;

hold off;

if debug; disp(['...Leaving ',mfilename]); end


