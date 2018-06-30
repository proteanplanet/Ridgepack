function [faxis,psd,conf,cutconf,description,units,V,maxlag,spec]=ridgepack_psd(ts,df,method,ts2)

% ridgepack_psd - Calculate a probability spectral density or cross-spectra
%
% function [faxis,psd,conf,cutconf,description,units,V,maxlag,spec]=ridgepack_psd(ts,df,method,ts2)
%
% This function generates a power spectra for the time series ts. The
% function's default uses Welch's method to calculate a periodogram 
% with degrees of freedom based on the amount of averaging used 
% for hamming-windowed data.  The size of the windows is provided 
% in terms of the effective frequency resolution of the final 
% spectra, df.  The larger the value of df, the greater confidence
% one can have in the spectra, but the greater the smoothing
% and decrease in frequency resolution of the final spectra.
% There is no overlap between windows. A periodogram can be plotted
% instead of a spectrogram by setting the method to 'periodogram' (the
% default is 'welch').
%
% Alternatively, the covariance method can be used to approximate
% the power spectral density using the Parzen window. In this case,
% the maximum lag to be used for the calculation is calculated 
% based on the specified df value.  This method uses a biased
% covariance estimate (see ridgepack_cov for more information), and
% can be specified by setting method to 'covariance' (see below).
% The covariance method is automatically used if cross spectra
% are requested.
%
% Note that there is a facility for second time series to be specified, 
% and if this is provided then the cross spectra of ts and ts2 is 
% calculated instead of the power spectra of ts.  If the cross spectra 
% is specified, then this function checks that ts and ts2 are 
% synchronized.
%
% If no output is specified, the call simply plots the PSD using,
% ncplotpsd otherwise output is given and the PSD is not plotted.
% If the argument sp is set to true a spectrogram is plotted, and 
% the output argument faxis is the axis handle for the spectrogram.
%
% The time series may be complex in order to calculate rotary 
% spectra, or else a simple one-sided spectra is calculated for
% real time series.
%
% INPUT:
%
% ts  - time series for power spectra
% df  - effective resolution of spectra (window size) in cycles per day.
%       Typical values are 0.0625 (relatively course) - 0.25  (smoothed)
%       If this argument is eliminated, the default is 0.125.
% method - method used to calcuate power spectra density. If method
%       is left blank, then the Welch averaging method is used,
%       otherwise if method is set to 'covariance', then the
%       covariance method is used with a Parsen window. Note that
%       if cross spectra are calculated, the covariance method
%       is automatically selected. If a periodogram instead of a
%       spectragram is required, then method should be set to
%       'periodogram'.
% ts2 - optional second time series for calculating cross-spectra.  This 
%       can either be omitted or entered as an empty matrix if cross-spectra
%       are not required.
%
%
% OUTPUT:
%
% faxis       - frequency axis of psd (Nx1) when specified with psd
%               or else if specified alone (one output argument) this 
%               is the axis handle on the spectrogram.
% psd         - power spectral density in dB relative to 1 (Nx1)
% conf        - lower and upper confidence bounds on psd (Nx2) 
%               (first index of second dimension is the lower error bound, 
%               second index of second dimension is the upper error bound)
%               Note that conf must be mulitplied by psd to give the upper
%               and lower bounds, and this is automatically done if using
%               ncplotpsd.
% cutconf     - confidence limit placed on data (preset to 95%)
% description - cell array of text with a description of windowing procedure
% units       - PSD units
% V           - Degrees of freedom of spectral estimate
% maxlag      - The maximum lag in days if using the covariance method 
% spec	      - The complex spectral values of the psd (equivalent to the
%               windowed fft of the covariance function).  This is only
%               available for the 'covariance' method specification.
%
% Ridgepack Version 1.0
% Andrew Roberts, Naval Postgraduate School, March 2018 (afrobert@nps.edu)
% Naval Postgraduate School, Department of Oceanography

global debug;
if debug; disp(['Entering ',mfilename,'...']); end

if ~isa(ts,'timeseries')
	error('The first argument is not a timeseries')
elseif nargin==1 % set default window size
	df=0.125;
end

crosspect=0;
maxlag=[];
spec=[];

if nargin<3
	method='welch';
	sp=0;
elseif nargin>3
	crosspect=1;
	sp=0;
	method='covariance';
elseif strcmpi(method,'covariance')
	method='covariance';
	sp=0;
elseif strcmpi(method,'periodogram')
	method='welch';
	sp=1;
elseif strcmpi(method,'welch')
	method='welch';
	sp=0;
end

% remove mean before proceeding
ts=detrend(ts,'constant');

% remove linear trend before proceeding
ts=detrend(ts,'linear');

% do the same for other timeseries
if crosspect
 ts2=detrend(ts2,'constant');
 ts2=detrend(ts2,'linear');
end

% check units of timeseries are days
if ~strcmp('days',ts.TimeInfo.Units)
	error('Time is not in units of days')
end

% set sampling interval from cleaned time series
if isnan(ts.TimeInfo.Increment)
        deltat=mean(ts.Time(2:end)-ts.Time(1:end-1));
        if debug
          disp(['Samples not uniformly distributed in time, mean interval is: ',...
                 num2str(mean(deltat))]);
        end
        %error('Reset the ts.TimeInfo.Increment variable or resample the timeseries');
else
	deltat=ts.TimeInfo.Increment; 
end

% synchronize cross spectra time series
if crosspect
 [ts ts2]=synchronize(ts,ts2,'uniform','interval',deltat);
 ts.TimeInfo.UserData=['dt=',num2str(deltat)]; 
 ts2.TimeInfo.UserData=['dt=',num2str(deltat)]; 
elseif strcmp(method,'covariance')
 ts2=ts;
end

if strcmp(method,'welch')

 % calculations using the Welch averaging method

 % start length
 Nstart=length(ts.data);

 % determine size of window based on effective frequency resolution df
 nwind=fix(1/(df*deltat)) ; % number of samples in the window
 disp(['Window size is ',num2str(nwind),' samples']);

 % use maximum possible time series length
 N=nwind*fix(length(ts.data)/nwind) ; 

 % remove samples that cannot be used from the end of the record
 ts=delsample(ts,'Index',[N+1:1:length(ts.data)]);
 ts.TimeInfo.Increment=deltat; 

 % inform user of critical numbers up to here
 disp([num2str(N),' samples used from total of ',num2str(Nstart)])

 % calculate power in dB/frequency bin with no overlap (hamming filtered)
 % using the pwelch command to calculate a spectra on a time series. 
 % Calculate rotary (two sided) spectra for complex timeseries.
 % There is no overlap between adjacent windows.
 if isreal(ts.Data)
  psd=pwelch(ts.Data,hamming(nwind),0,N,[],'onesided');
 else
  psd=fftshift(pwelch(ts.Data,hamming(nwind),0,N,[],'twosided'));
 end

 % Calculate the equivalent degrees of freedom based on Priestley, 1981, "Spectral
 % Analysis and Time Series", Academic Press: 2.5164*N/(window half width)
 % for a Hamming Window or "Tukey-Hamming Window" for a one-sided
 % power spectral density estimate. Here, the width of the window is nwind. 
 % Note that for a two sided spectra, the degrees of freedom is exactly 
 % half the single-sided spectra.
 if isreal(ts.Data)
  V=2*2.5164*N/nwind;  
 else
  V=2.5164*N/nwind;  
 end

elseif strcmp(method,'covariance')

 % calculate using the covariance method

 % get length of time series
 N=length(ts.Data);

 % calculate maximum lag based on requested frequency resolution    
 maxlag=fix(1/(2*df*deltat))
 if maxlag>(N/2)-1
  disp('Maximum lag exceeds (N/2)-1. Resetting to (N/2)-1')
  maxlag=(N/2)-1;
 end

 % calculate cross-covariance    
 crosscov=xcov(ts.Data,ts2.Data,maxlag,'biased');
 nwind=length(crosscov);

 % Calculate Parzen windowed power spectral density with one or two sided spectra 
 % depending on whether or not at least one of the time series is complex, otherwise 
 % make spectra one sided.  This method is based on  Priestley, 1981, "Spectral
 % Analysis and Time Series", Academic Press. 
 if isreal(ts.Data) & isreal(ts2.Data)
  psd=abs(fft(parzenwin(nwind).*crosscov));
  psd=2*psd(1:maxlag+1);
 else
  psd=abs(fftshift(fft(parzenwin(nwind).*crosscov)));
 end
 spec=fftshift(fft(parzenwin(nwind).*crosscov));

 % Calculate the degrees of freedom (check this is correct)
 lagwind=zeros([nwind 1]);
 for j=-maxlag:maxlag;
  lagwind(j+maxlag+1)=(nwind-abs(j))/nwind;
 end
 if isreal(ts.Data) & isreal(ts2.Data)
  V=2*length(ts.Data)/sum((lagwind.*parzenwin(nwind)).^2);
 else
  V=length(ts.Data)/sum((lagwind.*parzenwin(nwind)).^2);
 end

 % set length of ts.Data for the purpose of calculating the frequency axis
 ts=delsample(ts,'Index',nwind+1:N);
 ts.TimeInfo.UserData=['dt=',num2str(deltat)]; 

else

 error('Incorrect method specification')

end

% set the confidence interval (in %)
cutconf=95;

% set the fraction of the points that may fall outside the preferred
% level of confidence
alpha2=(1-cutconf/100)/2; 

% Now calculate the upper and lower confidence interval as using chi-squared
% distribution following Emery and Thomson, 2004, "Data Analysis Methods
% in Physical Oceanography", Elsevier. Note that conf must be muliplied by
% psd to give the upper and lower limits (i.e. psd*conf).
conf=V./chi2inv([1-alpha2,alpha2],V);

% generate frequency axis based on length of fft (N) using timeaxis.m;
[faxis,fs,df]=ridgepack_faxis(ts); faxis=faxis(end-length(psd)+1:end);

if nargout<7 & debug
 % inform the user of the frequency resolutions
 disp(['The Nyquist frequency is ',num2str(max(faxis)),' cycles per day'])
 disp(['The actual bandwidth is ',num2str(df),' cycles per day'])
end

% calculate effective frequency resolution based on length of windows
edf=1/(deltat*nwind);

% create description of the windowing process
if strcmp(method,'welch')
 description={['Effective PSD bandwidth ',num2str(edf,'%1.15'),' cycles {day$^{-1}$}'],...
              ['N = ',num2str(N),', ',num2str(deltat*24,3),' hour samples (samples/window = ',...
 	      num2str(nwind),')']};
elseif strcmp(method,'covariance')
 description={['Effective PSD bandwidth ',num2str(edf,'%1.5f'),' cycles {day$^{-1}$}'],...
              ['N = ',num2str(N),', ',num2str(deltat*24,3),' hour samples (max lag ',...
 	      num2str(maxlag*deltat,'%5.2f'),' days)']};
else
 error('Incorrect method specification')
end

% calculate units
if isempty(ts.DataInfo.Unit); units=''; else; units=['{(',ridgepack_units(ts),')$^{2}$}']; end

% If no output arguments are specified, plot the spectrum
if nargout==0 & ~sp 
	ridgepack_plotpsd(faxis,psd,{},conf,cutconf,description,units);
	title({ts.Name,ts.DataInfo.UserData},'HorizontalAlignment','center');
	clear faxis;
elseif nargout==0 & sp
	if isreal(ts.Data)
 	 S=ridgepack_db(abs(spectrogram(ts.Data,hamming(nwind),0,N,[])))';
	else
 	 S=ridgepack_db(abs(fftshift(spectrogram(ts.Data,hamming(nwind),0,N,[]),1)))';
        end  
	if min(faxis)<0; % wrap axes
         disp('Wrapping axes');
	 faxis=cat(1,-faxis(end),faxis);
	 S=cat(1,S(end,:),S(:,:));
	end
	imagesc(faxis,ts.Time(1:nwind:end),S);
	datetick('y','keeplimits');
	set(gca,'Layer','top');
	box on; 
	axis tight; 
	xlim([-4.5 4.5]);
	colorbar;
	colormap jet;
	xlabel(['Frequency (cycles {day$^{-1}$})']);
	title(ridgepack_cellcat({ts.Name,'(dB {bandwidth$^{-1}$})'}));
	clear faxis;
end

description

if debug; disp(['...Leaving ',mfilename]); end

