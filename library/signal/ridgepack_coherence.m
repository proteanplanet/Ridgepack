function [faxis,coh,cohou,conf,cutconf,description]=ridgepack_coherence(ts1,ts2,df,method)

% ridgepack_coherence - Generates spectral coherence squared for time series
%
% function [faxis,coh,cohou,conf,cutconf,description]=ridgepack_coherence(ts1,ts2,df,method)
%
% This function generates spectral coherence squared for time series. This function
% uses Welch's method to calculate coherence using non-overlapping hamming-windowed
% data. There are three potential types of spectral coherence that can be generated:
%
% 1) Standard spectral coherence for two real time series, ts1 and ts2
%
% 2) Rotary spectral coherence for two complex time series, ts1 and ts2
%    where the real and imaginary part of each represent orthogonal
%    components of two dimensional vector. The output from this provide
%    both the inner (co-rotating) coherence and outer (contra-rotating)
%    coherence.
%
% 3) Auto coherence of a single, complex time series ts1 with real and imaginary
%    components representing the orthogonal components of a 2D vector.
%    The output indicates the degree to which there is spectral coherence
%    between the oppositely rotating components of the vector.
%
% If no output arguments are provided, this function plots the coherence 
% with the ncplotcoh.m function. Note that phase information is not provided
% as part of the cross-spectral coherence because it is typically difficult
% to obtain useful information from this parameter.
%
% INPUT:
%
% ts1    - Matlab time series object (real or complex)
% ts2    - Matlab time series object (real or complex) spanning an overlapping
%          period of time with ts1. For auto coherence, this can be entered
%          as an empty matrix.
% df     - Effective resolution of spectra (window size) in cycles per day.
%          Typical values are 0.0625 (relatively course) - 0.25  (smoothed)
%          If this argument is eliminated, the default is 0.125.
% method - Defines the method used to calculate coherence.  The default
%          method is the covariance method using a biased covariance
%          estimator and Parzen window.  By specifying the Welch method
%          the default can be changed to use window averaging. Hence method 
%          is either omitted for the default, or set to 'welch' or 'covariance'.
%
%
% OUTPUT:
%
% faxis       - frequency axis of psdb (Nx1)
% coh         - inner coherence, standard spectral coherence or auto coherence (Nx1)
% cohou       - outer coherence
% conf        - confidence array for each coherence value (Nx1)
% cutconf     - confidence limit placed on data (preset to 95%)
% description - cell array of text with a description of windowing procedure
%
% Ridgepack Version 1.0
% Andrew Roberts, Naval Postgraduate School, March 2018 (afrobert@nps.edu)
% Naval Postgraduate School, Department of Oceanography
%

% check input 
if nargin<2
	error('There must be three input arguements')
elseif nargin<3; 
	df=0.125; 
end

if nargin<4
	method='covariance';
elseif strcmpi(method,'welch')
	method='welch';
elseif strcmpi(method,'covariance') 
	method='covariance';
else
	error('incorrect method specification')
end

if ~isa(ts1,'timeseries')
	error('ts1 must be a timeseries')
elseif ~isa(ts2,'timeseries')
	auto=true;
	disp('Calculating Autocoherence')
	ts2=ts1;
	ts2.Data=conj(ts1.Data);
        ts2.Timeinfo.UserData=ts1.Timeinfo.UserData;
else
	auto=false;
end

% remove mean before proceeding
ts1=detrend(ts1,'constant'); ts2=detrend(ts2,'constant');

% remove linear trend
ts1=detrend(ts1,'linear'); ts2=detrend(ts2,'linear');

% check units of timeseries are days
if ~strcmp('days',ts1.TimeInfo.Units) | ~strcmp('days',ts2.TimeInfo.Units)
 error('Time is not in units of days for one or both of the time series')
end

% set sampling interval from cleaned time series
if ~isempty(ts1.TimeInfo.Increment) & ischar(ts1.TimeInfo.UserData) & ...
    length(ts1.TimeInfo.UserData)>3 & strcmpi(ts1.Timeinfo.UserData(1:3),'dt=')
 try
  deltat=str2num(ts1.TimeInfo.UserData(4:end));
 catch
  error('Unable to convert ts.Timeinfo.UserData to a time interval')
 end
elseif isnan(ts1.TimeInfo.Increment)  
 deltat=mean(ts1.Time(2:end)-ts1.Time(1:end-1));
 disp(['Samples not uniformly distributed in time for ts1, mean interval is: ',num2str(mean(deltat))]);
 %error('Reset the ts1.TimeInfo.Increment variable or resample the timeseries');
elseif isnan(ts2.TimeInfo.Increment) 
 deltat=mean(ts2.Time(2:end)-ts2.Time(1:end-1));
 disp(['Samples not uniformly distributed in time for ts2, mean interval is: ',num2str(mean(deltat))]);
 %error('Reset the ts2.TimeInfo.Increment variable or resample the timeseries');
else
 deltat=ts1.TimeInfo.Increment;
end

% synchronize datasets to the first timeseries
[ts1 ts2]=synchronize(ts1,ts2,'uniform','interval',deltat);
ts1.Timeinfo.UserData=['dt=',num2str(deltat)];
ts2.Timeinfo.UserData=['dt=',num2str(deltat)];

% set variables specific to the method being used
if strcmp(method,'welch')

 % start length
 Nstart=length(ts1.Data);

 % number of samples in each window
 nwind=fix(1/(df*deltat)) ; 
 disp(['Window size is ',num2str(nwind),' samples']);

 % use maximum possible time series length
 N=nwind*fix(length(ts1.Data)/nwind) ; 

 % remove samples that cannot be used from the end of the record
 ts1=delsample(ts1,'Index',[N+1:1:length(ts1.Data)]);
 ts2=delsample(ts2,'Index',[N+1:1:length(ts2.Data)]);

 % display information
 disp([num2str(N),' samples used from total of ',num2str(Nstart)])

elseif strcmp(method,'covariance')

 % get length of time series
 N=length(ts1.Data);

 % calculate maximum lag based on requested frequency resolution    
 maxlag=fix(1/(2*df*deltat));
 if maxlag>(N/2)-1
  disp('Maximum lag exceeds (N/2)-1. Resetting to (N/2)-1')
  maxlag=(N/2)-1;
 end
 nwind=(2*maxlag)+1;

 % set N to actual nwind length
 N=nwind;

else

 error('Incorrect method specification')

end

% set sample rates explicitly
ts1.Timeinfo.UserData=['dt=',num2str(deltat)];
ts2.Timeinfo.UserData=['dt=',num2str(deltat)];


% calculate power in dB/frequency bin with 0 overlap (hamming)
if isreal(ts1.Data) & isreal(ts1.Data)

	disp('Calculating scalar coherence');

	if strcmp(method,'welch')

	 % scalar magnitude squared coherence with no window overlap
	 coh=mscohere(ts1.Data,ts2.Data,hamming(nwind),0,N,'onesided');

	 % Calculate the equivalent degrees of freedom based on Priestley, 1981, "Spectral
	 % Analysis and Time Series", Academic Press: 2.5164*N/(window half width)
	 % for a Hamming Window or "Tukey-Hamming Window" for a one-sided
	 % power spectral density estimate. Here, the width of the window is nwind. 
	 V=2*2.5164*N/nwind;

	elseif strcmp(method,'covariance')

         % calculate cross spectra from covariance power spectral density 
	 [faxis,crossspectra,conf,cutconf,description,units,V,maxlag]=ridgepack_psd(ts1,df,'covariance',ts2);
	 [faxis,ts1spectra,conf,cutconf,description,units,V,maxlag]=ridgepack_psd(ts1,df,'covariance');
	 [faxis,ts2spectra,conf,cutconf,description,units,V,maxlag]=ridgepack_psd(ts2,df,'covariance');

	 coh=(crossspectra.^2)./(ts1spectra.*ts2spectra);

        else

         error('Incorrect method specification')

        end


	% no counter-rotating components
	cohou=[];

elseif auto

	disp('Calculating autocoherence');

	if strcmp(method,'welch')

	 % coh for the autocoherence case the coherence for counted-rotating components 
	 % of the same timeseries, or the outer-coherence of a timeseries with itself.
	 coh=mscohere(ts1.Data,conj(ts2.Data),hamming(nwind),0,N,'twosided');

	 % Calculate the degrees of freedom. In this case, becuase both the positive
         % and negative frequencies have the same coherence, the degrees of freedom
         % doubles for the positive side of the coherence being plotted.	
	 V=2*2.5164*N/nwind;

	elseif strcmp(method,'covariance')

         % calculate cross spectra from covariance power spectral density 
	 [faxis,crossspectra,conf,cutconf,description,units,V,maxlag]=ridgepack_psd(ts1,df,'covariance',ts2);
	 [faxis,ts1spectra,conf,cutconf,description,units,V,maxlag]=ridgepack_psd(ts1,df,'covariance');
	 [faxis,ts2spectra,conf,cutconf,description,units,V,maxlag]=ridgepack_psd(ts2,df,'covariance');
	 coh=(crossspectra.^2)./(ts1spectra.*ts2spectra);

        else

         error('Incorrect method specification')

        end

	% Both sides of an autocoherence spectra are the same, so there is only a need
	% to keep the positive side of the cross-coherence. 
	if floor(N/2)==N/2 % even N
	 coh=coh(1:N/2+1);
	else % odd N
	 coh=coh(1:(N+1)/2);
	end

	% no counter-rotating component
	cohou=[];

elseif ~isreal(ts1.Data) | ~isreal(ts1.Data)

	disp('Calculating rotary coherence');

	if strcmp(method,'welch')

	 % coh is the coherence for co-rotating components with no window overlap
	 coh=fftshift(mscohere(ts1.Data,ts2.Data,hamming(nwind),0,N,'twosided'));

	 % cohou is the coherence for counter-rotating components (the conjugate
	 % of the second timeseries sets the rotation in the opposite direction)
	 cohou=fftshift(mscohere(ts1.Data,conj(ts2.Data),hamming(nwind),0,N,'twosided'));

	 % Calculate the degrees of freedom (see above for explanation). Note that for 
	 % a two sided spectra, the degrees of freedom is exactly half the single-sided spectra.
	 V=2.5164*N/nwind;

	elseif strcmp(method,'covariance')

         % calculate cross spectra from covariance power spectral density 
	 [faxis,crossspectra,conf,cutconf,description,units,V,maxlag]=ridgepack_psd(ts1,df,'covariance',ts2);
	 [faxis,ts1spectra,conf,cutconf,description,units,V,maxlag]=ridgepack_psd(ts1,df,'covariance');
	 [faxis,ts2spectra,conf,cutconf,description,units,V,maxlag]=ridgepack_psd(ts2,df,'covariance');
	 coh=(crossspectra.^2)./(ts1spectra.*ts2spectra);

         % calculate cross spectra from covariance power spectral density 
	 ts2.Data=conj(ts2.Data);
	 [faxis,crossspectra,conf,cutconf,description,units,V,maxlag]=ridgepack_psd(ts1,df,'covariance',ts2);
	 [faxis,ts1spectra,conf,cutconf,description,units,V,maxlag]=ridgepack_psd(ts1,df,'covariance');
	 [faxis,ts2spectra,conf,cutconf,description,units,V,maxlag]=ridgepack_psd(ts2,df,'covariance');
	 cohou=(crossspectra.^2)./(ts1spectra.*ts2spectra);
	 ts2.Data=conj(ts2.Data);

        else

         error('Incorrect method specification')

        end

else

	error('One timeseries is a vector (complex) series, the other is scalar')

end

% calculate confidence interval (V=equivalent degrees of freedom)
% This is based on equations provided in Emery and Thomson, 2004, 
% "Data Analysis Methods in Physical Oceanography", Elsevier.
cutconf=99.0; % 99% confidence interval
conf=ones(size(coh))*(1-((1-cutconf/100)^(1/(V-1))));

% calculate effective frequency resolution based on length of (nwind);
edf=1/(deltat*nwind);

% create description of the windowing process
if strcmp(method,'welch')

 if auto
  description={['Effective bandwidth = ',num2str(edf),' cycles {day$^{-1}$}'],...
               ['N = ',num2str(N),', ',num2str(deltat*24,3),' hour samples (window = ',...
  	       num2str(nwind),')'],'Autocoherence'};
 else
  description={['Effective bandwidth = ',num2str(edf),' cycles {day$^{-1}$}'],...
               ['N = ',num2str(N),', ',num2str(deltat*24,3),' hour samples (window = ',...
  	       num2str(nwind),')']};
 end
elseif strcmp(method,'covariance')

 % shorten ts1 for calculation of correct df using faxis
 ts1=delsample(ts1,'Index',nwind+1:N); 
 ts1.Timeinfo.UserData=['dt=',num2str(deltat)];

 if auto
  description={['Bandwidth = ',num2str(edf),' cycles {day$^{-1}$}'],...
               ['N = ',num2str(N),', ',num2str(deltat*24,3),' hour samples (maximum lag = ',...
  	       num2str(maxlag*deltat),' days)'],'Autocoherence'};
 else
  description={['Bandwidth = ',num2str(edf),' cycles {day$^{-1}$}'],...
               ['N = ',num2str(N),', ',num2str(deltat*24,3),' hour samples (maximum lag = ',...
  	       num2str(maxlag*deltat),' days)']};
 end
else
 error('Incorrect method specification')
end

% CHANGE ANDREW ROBERTS JUNE 18 2011 - ERROR IN THIS SECTION WHICH IS REDUNDENT ANYWAY
% generate frequency axis based on length of fft (N) using timeaxis.m;
%[faxis,fs,df]=ridgepack_faxis(ts1); 
%faxis=faxis(end-length(coh)+1:end);

% inform the user of the frequency resolutions
disp(['The Nyquist frequency is ',num2str(max(faxis)),' cycles per day'])
disp(['The actual bandwidth is ',num2str(max(df)),' cycles per day'])

% if no output arguments are specified, plot the spectrum
if nargout==0
	if isreal(ts1.Data) | isreal(ts2.Data)
	  ncplotcoh(faxis,coh,cohou,{},conf,cutconf,description,0,min(max(faxis)));
        else
	  ncplotcoh(faxis,coh,cohou,{},conf,cutconf,description,max(min(faxis),-4.5),min(max(faxis),4.5));
        end
	if auto ; 
		title(ts1.Name); 
	elseif isreal(ts1.Data);
		title({ts1.Name,ts2.Name}); 
	end
	clear faxis;
end

end

