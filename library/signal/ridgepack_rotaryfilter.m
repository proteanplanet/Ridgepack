function [Wns_clock,Wns_anticlock,ts_filt,scale,period,gain,faxis,pband]=...
            ridgepack_rotaryfilter(ts,band,m)

% ridgepack_rotaryfilter - performs a rotary wavelet transform on complex timeseries
%
% function [Wns_clock,Wns_anticlock,ts_filt,scale,period,gain,faxis,pband]=...
%             ridgepack_rotaryfilter(ts,band,m)
%
% This function performs a rotary wavelet transform on complex timeseries, and provide 
% reconstructions of partial components of timeseries, using a 18th order DOG wavelet by 
% default.
%
% Inputs:
% ts       - matlab timeseries object
% passband - Two element vector giving the lower and upper pass frequencies in cycles/day
% m        - filter order (an integer number betweeen 6 and 50 where m-2 is divisible by 4) 
%            {optional}
%
% Outputs:
% Wns_clock     - clockwise rotary power for the tolerance stipulated
% Wns_anticlock - anticlockwise rotary power for the tolerance stipulated
% ts_filt       - matlab timeseries collection of reconstructed and filtered time series
% scale         - wavelet scale for the above Wns outputs
% period        - associated period for each scale
% gain          - filter gain in dB as a structure for each timeseries in ts_filt
% faxis         - frequency axis in cycles/day of the gain signal
% pband         - transition band required to meet filter tolerance in cycles/day
%
% Ridgepack Version 1.0
% Andrew Roberts, Naval Postgraduate School, March 2018 (afrobert@nps.edu)
% Naval Postgraduate School, Department of Oceanography

% copy timeseries
tscopy=ts;

% get the timeseries timestep
[ts,dt]=ridgepack_uniform_sample_check(ts);

% prepare time series by removing DC signal and linear trend
if nargin>0
 z1=ts.data;
 ts=detrend(ts,'constant');
 ts=detrend(ts,'linear');
 meanval=z1-ts.data;
 z1=ts.data;
else
 error('No inputs')
end

% sort the band correctly if needed;
if nargin==1
 band=[0 24];
else
 if ~(length(band)==2) 
  error('Band length is incorrect');
 elseif any(band<0)
  error('Band frequencies must be greater than zero');
 end
 band=sort(band);
 if diff(band)==0
  error('band must be two different frequencies')
 end
end

% Check to see whether meanvalue should be added to the Total timeseries.
% (i.e. if the clockwise+anticlockwise band limited signal extends from freq=0)
if any(band==0)
 addmean=true;
else
 addmean=false;
end

% Diagnostics
bandwidthtext=[num2str(band(1)),' to ',num2str(band(2)),' cycles/day'];
disp(['Filtering from ',bandwidthtext]); 

% get timeseries length
N=length(ts.data);

% get fourier transform of timeseries
xk=fft(z1);

% construct angular frequency array
n=[0:1:N-1]';
wk=[(1-N)/2:1:(N-1)/2]*2*pi/(N*dt);
midpoint=fix(length(wk)/2);
wk=[wk(midpoint+1:end) wk(1:midpoint)];

% set frequency resolution
df=diff(wk(1:2));

% filter tolerance (nominally 1dB) and spectral resolution for iterations 
% and -10dB for transition area.
%Tolerance=1;
Tolerance=2;


% set DOG wavelet order
if nargin<3
 m=18;
 disp(['Filter order preset to m=',num2str(m)])
elseif m<6
 m=6;
 disp(['Changing filter order to m=',num2str(m)])
end

% construct scale array
dj=0.025; % scale increment
s0=dt*sqrt(m+0.5)/pi;
maxscale=fix(log2(N*dt/s0)/dj);
scale=s0*2.^((0:maxscale)*dj);


% check that the filter order produces a positive real Fourier Transform and does
% not exceed computational limits (i.e. the width of the DOG is within range).
diagnostic_output=0;

if floor((m-2)/4)~=(m-2)/4
 m=4*floor((m-2)/4)+2;
 s0=dt*sqrt(m+0.5)/pi;
 maxscale=fix(log2(N*dt/s0)/dj);
 scale=s0*2.^((0:maxscale)*dj);
 diagnostic_output=1;
end

maxcomp=((max(scale)*max(wk)).^m).*exp(-((max(scale)*max(wk)).^2)/2)./sqrt(gamma(m+0.5));
while isnan(maxcomp)
 m=m-4;
 if m<6; error('CODE ERROR: m<6'); end
 s0=dt*sqrt(m+0.5)/pi;
 maxscale=fix(log2(N*dt/s0)/dj);
 scale=s0*2.^((0:maxscale)*dj);
 maxcomp=((max(scale)*max(wk)).^m).*exp(-((max(scale)*max(wk)).^2)/2)./sqrt(gamma(m+0.5));
 diagnostic_output=1;
end

if diagnostic_output; disp(['Reducing filter order to m=',num2str(m)]); end

% set maximum iterations
maxalpha=100;

% initialize wavelet parameters and matrices
period=zeros(size(scale)); 
Wns=zeros(size(scale'*wk));
Wns_clockwise=zeros(size(Wns));
Wns_anticlockwise=zeros(size(Wns));

% plot demonstration fourier transform of wavelet corresponding to 
% Torrence and Compo (1998) right panel Figure 2(d)
demo=0;
if demo 

 clf

 scale_demo=10*dt;
 psi_hat_clockwise=-(sqrt(2*pi*scale_demo/dt))*((j.^m)/sqrt(gamma(m+0.5))).*...
                    (((scale_demo*(wk.*(wk<0))).^m).*exp(-((scale_demo*wk).^2)/2));
 psi_hat_anticlockwise=-(sqrt(2*pi*scale_demo/dt))*((j.^m)/sqrt(gamma(m+0.5))).*...
                    (((scale_demo*(wk.*(wk>0))).^m).*exp(-((scale_demo*wk).^2)/2));
 abscissa=[(1-N)/2:1:(N-1)/2]*scale_demo/(N*dt);
 clf
 h(1)=plot(abscissa,psi_hat_clockwise([midpoint+1:end 1:midpoint]),'b')
 hold on
 h(4)=plot(abscissa,psi_hat_anticlockwise([midpoint+1:end 1:midpoint]),'r')

 MM=m;
 m=6;
 psi_hat_clockwise=-(sqrt(2*pi*scale_demo/dt))*((j.^m)/sqrt(gamma(m+0.5))).*...
                    (((scale_demo*(wk.*(wk<0))).^m).*exp(-((scale_demo*wk).^2)/2));
 psi_hat_anticlockwise=-(sqrt(2*pi*scale_demo/dt))*((j.^m)/sqrt(gamma(m+0.5))).*...
                    (((scale_demo*(wk.*(wk>0))).^m).*exp(-((scale_demo*wk).^2)/2));
 abscissa=[(1-N)/2:1:(N-1)/2]*scale_demo/(N*dt);
 h(2)=plot(abscissa,psi_hat_clockwise([midpoint+1:end 1:midpoint]),'b--')
 h(5)=plot(abscissa,psi_hat_anticlockwise([midpoint+1:end 1:midpoint]),'r--')

 m=30;
 psi_hat_clockwise=-(sqrt(2*pi*scale_demo/dt))*((j.^m)/sqrt(gamma(m+0.5))).*...
                    (((scale_demo*(wk.*(wk<0))).^m).*exp(-((scale_demo*wk).^2)/2));
 psi_hat_anticlockwise=-(sqrt(2*pi*scale_demo/dt))*((j.^m)/sqrt(gamma(m+0.5))).*...
                    (((scale_demo*(wk.*(wk>0))).^m).*exp(-((scale_demo*wk).^2)/2));
 abscissa=[(1-N)/2:1:(N-1)/2]*scale_demo/(N*dt);
 h(3)=plot(abscissa,psi_hat_clockwise([midpoint+1:end 1:midpoint]),'b:')
 h(6)=plot(abscissa,psi_hat_anticlockwise([midpoint+1:end 1:midpoint]),'r:')

 plot([0 0],[0 6],'k--')
 plot([0 0],[0 6],'k--')
 xlim([-2 2])
 xlabel('s \omega / (2 \pi)')
 set(gca,'box','off')

 legend(h(4:6),{['anticlockwise m = ',num2str(MM)],['anticlockwise m = 6'],['anticlockwise m = 30']})
 legend('boxoff')
 set(gca,'box','on')

 pos=get(gca,'position');
 ah=axes('position',pos,'visible','off');
 legend(ah,h(1:3),{['clockwise m = ',num2str(MM)],['clockwise m = 6'],['clockwise m = 30']},'location','northwest')

 legend('boxoff')

 ncfprint('png','Demonstration_psi_hat',1,1)
 ncfprint('epsc','Demonstration_psi_hat',1,1)
 hold off

 m=MM;

 return
end


% do wavelet transform using DOG mother wavelet for each scale
for i=1:length(scale)

 % rotary fourier transform of DOG wavelet 
 psi_hat_clockwise=-(sqrt(2*pi*scale(i)/dt))*((j.^m)/sqrt(gamma(m+0.5))).*(((scale(i)*(wk.*(wk<0))).^m).*exp(-((scale(i)*(wk)).^2)/2));

 psi_hat_anticlockwise=-(sqrt(2*pi*scale(i)/dt))*((j.^m)/sqrt(gamma(m+0.5))).*(((scale(i)*(wk.*(wk>0))).^m).*exp(-((scale(i)*(wk)).^2)/2));

 psi_hat=-(sqrt(2*pi*scale(i)/dt))*((j.^m)/sqrt(gamma(m+0.5))).*(((scale(i)*wk).^m).*exp(-((scale(i)*wk).^2)/2));

 % fourier wavelength
 period(i)=2*pi*scale(i)/sqrt(m+0.5);

 % clockwise, anticlockwise and full wavelet transform
 Wns_clockwise(i,:)=ifft(xk.*conj(psi_hat_clockwise)');
 Wns_anticlockwise(i,:)=ifft(xk.*conj(psi_hat_anticlockwise)');
 Wns(i,:)=ifft(xk.*conj(psi_hat)');

end

% convert bands to periods and keep within represented periods
band=1./band;
band(1)=max(min(band(1),max(period(:))),min(period(:)));
band(2)=min(max(band(2),min(period(:))),max(period(:)));

% calculate complete band limited reconstruction
s=1./sqrt(scale);

% set PSD method
method='covariance';

% Determine empirical Cdelta*Psi0 function using Parseval's Theorem
[faxis,psd_input]=ridgepack_psd(ts,df,method);
ts.data=[dj*sqrt(dt)*(s*conj(Wns))]'+meanval;
[faxis,psd_output]=ridgepack_psd(ts,df,method);
p1=sum(psd_input);
p2=sum(psd_output);
CdeltaPsi0=sqrt(p2/p1);
disp(['Broadband CdeltaPsi0 = ',num2str(CdeltaPsi0)])

% calculate fully reconstructed field
if all(abs(z1)==0)
 ts.data=z1+meanval;
else
 ts.data=[dj*sqrt(dt)*(s*conj(Wns))/CdeltaPsi0]'+meanval;
end
ts.name='Reconstructed';
ts.DataInfo.UserData='Complete reconstructed signal with mean';
ts_filt.Name=['Complete broadband signal reconstruction with mean'];
ts_filt=tscollection(ts);

% calculate fully reconstructed gain
[faxis,psd]=ridgepack_psd(ts,df,method);
gain.Reconstructed=10*log10(psd./psd_input);

% determine periods that fit the bands
periodpass=find(period>=band(2) & period<=band(1));
freqpass=find(faxis>=1./band(1) & faxis<=1./band(2));

% calculate band-limited scale
s=1./sqrt(scale(periodpass));

% Determine empirical Cdelta*Psi0 function using Parseval's Theorem (band limited)
Wns=Wns(periodpass,:);
if all(abs(z1)==0)
 ts.data=z1;
else
 ts.data=[dj*sqrt(dt)*(s*conj(Wns))]';
end
[faxis,psd]=ridgepack_psd(ts,df,method);
p1=sum(psd_input(freqpass));
p2=sum(psd(freqpass));
CdeltaPsi0=sqrt(p2/p1);
disp(['Band limited CdeltaPsi0 = ',num2str(CdeltaPsi0)])

% calculate complete band limited reconstruction, adding the mean
% value if the filter bands include 0 cycles/day.
if addmean & all(abs(z1)==0)
 ts.data=meanval;
elseif all(abs(z1)==0)
 ts.data=z1;
elseif addmean
 ts.data=[dj*sqrt(dt)*(s*conj(Wns))/CdeltaPsi0]'+meanval;
else
 ts.data=[dj*sqrt(dt)*(s*conj(Wns))/CdeltaPsi0]';
end
ts.name='Total';
ts.DataInfo.UserData='Complete band-limited signal';
ts_filt.Name=['Complete band-limited signal reconstruction: ',bandwidthtext];
ts_filt=addts(ts_filt,ts);

% calculate reconstructed field gainf for combined CW and ACW signal
[faxis,psd]=ridgepack_psd(ts,df,method);
gain.Total=10*log10(psd./psd_input);

% calculate anticlockwise reconstructed signal
Wns_anticlock=Wns_anticlockwise(period>=band(2) & period<=band(1),:);
if all(abs(z1)==0)
 ts.data=z1;
else
 ts.data=[dj*sqrt(dt)*(s*conj(Wns_anticlock))/CdeltaPsi0]';
end
ts.name='Anticlockwise';
ts.DataInfo.UserData=['Band limited anticlockwise signal: ',bandwidthtext];
ts_filt=addts(ts_filt,ts);

% calculate anticlockwise reconstructed gain
[faxis,psd]=ridgepack_psd(ts_filt.Anticlockwise,df,method);
gain.Anticlockwise=10*log10(psd./psd_input);

% calculate clockwise reconstructed signal
Wns_clock=Wns_clockwise(period>=band(2) & period<=band(1),:);
if all(abs(z1)==0)
 ts.data=z1;
else
 ts.data=[dj*sqrt(dt)*(s*conj(Wns_clock))/CdeltaPsi0]';
end
ts.name='Clockwise';
ts.DataInfo.UserData=['Band limited clockwise signal: ',bandwidthtext];
ts_filt=addts(ts_filt,ts);

% calculate clockwise reconstructed gain
[faxis,psd]=ridgepack_psd(ts,df,method);
gain.Clockwise=10*log10(psd./psd_input);

disp(['Pass tolerance: ',num2str(Tolerance),' dB'])
disp(['  Maximum gain: ',num2str(max(gain.Clockwise)),' dB'])
disp(['  Minimum gain: ',num2str(min(gain.Clockwise)),' dB'])

% tolerance test for non-zero timeseries 
if all(abs(z1)==0)
 pband(1)=0;
 pband(2)=0;
elseif isempty(find(abs(gain.Clockwise)<Tolerance))
 clf
 plot(faxis,gain.Clockwise)
 title('Poorly configured filter')
 xlabel('Frequency')
 ylabel('Gain')
 error('Poorly configured filter')
else
 pband(1)=min(abs(faxis(abs(gain.Clockwise)<Tolerance)));
 pband(2)=max(abs(faxis(abs(gain.Clockwise)<Tolerance)));
end

disp(['     Pass band: ',num2str(pband(1)),' to ',num2str(pband(2))])

if max(gain.Clockwise)>Tolerance
 disp(['WARNING: Filter exceeds Tolerance of +/-',num2str(Tolerance),' dB'])
end

% add tolerance to gain
gain.Tolerance=Tolerance*ones(size(faxis));

% Output diagnostics about full reconstruction
disp('-------------------------------------------')
disp(['Filter Order: m=',num2str(m)]); 

diffs=tscopy.Data-ts_filt.Reconstructed.Data;
N=size(diffs);

RMS1=sqrt((1/N(1))*sum(tscopy.Data.^2));
RMS2=sqrt((1/N(1))*sum(ts_filt.Reconstructed.Data.^2));
RMSE=sqrt((1/N(1))*sum(diffs.^2));

disp(['RMS Original: ',num2str(abs(RMS1)),' m/s'])
disp(['RMS Reconstructed: ',num2str(abs(RMS2)),' m/s'])
disp(['RMSE: ',num2str(abs(RMSE)),' m/s'])
disp(['RMSE: ',num2str(100*abs(RMSE)/abs(RMS1)),'%'])

var1=var(tscopy.Data);
var2=var(ts_filt.Reconstructed.Data);
disp(['Variance Original: ',num2str(var1)])
disp(['Variance Reconstructed: ',num2str(var2)])
disp(['Variance Retained: ',num2str(100*var2/var1),'%'])

disp('-------------------------------------------')


