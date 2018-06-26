function [ts,b,a,description]=ridgepack_filter1(ts,Wn,kind,type,n,bandwidth,smooth)

% ridgepack_filter1 - FIR and IIR time series filtering routine 
%
% function [ts,b,a,description]=ridgepack_filter1(ts,Wn,kind,type,n,bandwidth,smooth)
% 
% This  function filters data along the time dimension of a time series
% object, ts, using a hamming-windowed linear-phase FIR digital filter or
% an IIR Butterworth filter.  
%
% The data in ts may be multi dimensional, but only the time axis of the data 
% is filtered. The filter order, n, can be specified, or, if n is left empty
% this function finds the minimum order,n, meeting the design criteria set by 
% the filter bandwidth, which is the transition band between passed and stopped 
% frequencies. For FIR filters, the stop band attenuation is preset at 30dB 
% (0.1% leakage or less in stop bands) and the pass band ripple is set at 99% or 
% more of the input signal. If the bandwidth is omitted, the transition bands 
% are kept within 0.5 cycles/day. For IIR filters, the stop band attenuation is 
% preset at 20dB (1% leakage or less in stop bands) and the pass band ripple 
% is set at 99% or more of the input signal. If the bandwidth is omitted, the 
% transition bands are kept within 1.0 cycles/day.
%
% The cutoff frequencies are specified in Wn, and must be specified in cycles/day.
%
% If no output arguments are specified, this function acts as a tool for filter
% design, and plots the frequency and impulse response functions using the 
% ncfreqzplot. Otherwise the output is a filtered time series, with the option 
% of also obtaining the difference equation coefficients b and a for the filter. 
% These coefficients can the be used in other programs as a pre-designed filter.  
%
% Although the filter design is causal, the filter is implemented as a 
% non-causal filter (centered on the data point being output) so that there 
% is zero phase shift (no delay) in the output. This means that the final set 
% of samples, at time t, incorporate values at time t+deltat and t-deltat rather 
% than just t-2*detlat as is the case for a causal filter with 2*deltat delay. 
% As such, the final filter applied is in fact non-causal.
%
% Because the filter has been carefully designed to distinguish pass, transition 
% and stop bands, there is no recoloring applied to the pass bands.  However
% this may be required if on chooses to ignore the present passband ripple
% criteria, because this could result in the power of the passband being 
% significantly increased or decreased above or below the input signal.
%
% INPUT:
%
% ts        - time series object generated using time series
%             or ridgepack_2timeseries within this package.
%
% Wn        - Cutoff frequency w1, or a vector giving the stop or pass
%             band filters [w1 w2 w3 ... wn]. The units of the frequency
%             must be in cycles per day.  If using the FIRH filter kind,
%             Wn specifies the edge of the passband (see fir1hibler for
%             more details).
% 
% kind      - This chooses the type of filter to be used:
%             FIR  - Classic Finite Impulse Response filter
%             FIRH - Cosine Finite Impulse Response filter (low pass only)
%             IIR  - Butterworth Infinite Impulse Response 
%  
%             The cosine FIRH filter uses the technique in:
%             Hibler, W. D., III, 1972: 
%             Removal of Aircraft motion from LASER profilometers, 
%             J. Geophys. Res., 77, 36.
%             
% type      - 'low','high','stop' or 'pass' pass filters.  The default
%             is a low pass filter if this argument is omitted. For
%             the 'stop' and 'pass' band filters, Wn should be a vector
%             [w1 w2 w3 w4 w5] with the first band 0<w<w1 being stopped
%             or passed, respectively, and the second band w1<w<w2 being
%             passed or stopped, respectively.
%
% n         - filter order. If this is left empty [] or omitted, the
%             function will find the optimal order for the given
%             specifications. If this value is provided, then optimal
%             n will not be found, but instead the requested 
%             filter will be provided regardless of whether or
%             not the design criteria are met.
%
% bandwidth - filter bandwidth (ws-wp) where ws is the stop
%             frequency and wp is the pass frequency for any
%             transition with cutoff wc=(wp+ws)/2 for any of 
%             wc=w1, w2, ..., wn. Default is 0.5 cycles/day if 
%             for FIR and FIRH filters, and 1 cycle/day for 
%             the IIR filter. This value may be omitted or 
%             entered as an empty [] value for the defaults.
%
% smooth    - smoothly fades in data for low and pass band filters
%             from a starting constant value specified as 'smooth'.
%	      This is a single number that is used as a constant starting 
%             field from which the variable field being filtered is
%             fased in using a sinuisoid with frequency w1/4. This 
%             option is useful to avoid impulsive initialization of models.
%            
%
%
% OUTPUT:
%
% ts        - Zero phase shift filtered Matlab time series object.
%
% b         - b coefficients of the convolution filter.
% 
% a         - a coefficient included for completeness (a==1 for 
%             all cases in this function).
%
% description - text giving a description of the filter.
%
% Ridgepack Version 1.0
% Andrew Roberts, Naval Postgraduate School, March 2018 (afrobert@nps.edu)
% Naval Postgraduate School, Department of Oceanography
%

% set initial choice values
firc=false;
firh=false;
iirb=false;

% check inputs
if nargin<2
	error('Must provide the time series and frequency cutoff');
elseif ~isa(ts,'timeseries')
	error('Input ts is not a time series');
elseif nargin<3
	kind='FIR';	
	type='low';
elseif nargin<4
	type='low';
end

if strcmp(type,'stop') & length(Wn)==1
	error('must specify at least two frequencies for stop band filter');
elseif strcmp(type,'stop') & length(Wn)==2
	type='stop';
elseif strcmp(type,'stop') 
	type='DC-0';
elseif strcmp(type,'pass') & length(Wn)<3
	error('Must specify at least three frequencies for pass band')
elseif strcmp(type,'pass')
	type='DC-1';
elseif ~strcmp(type,{'low','high'})
	error('type must be set to one of: low, high, stop, pass');
end

% get time and frequency axis information
[faxis,fs,df]=ridgepack_faxis(ts);

% set tolerances
disp('-------------------------------');
if nargin<5; n=[]; emptyn=true; elseif isempty(n); emptyn=true; else; emptyn=false; end

if strcmp(kind,'FIR') 

 disp('Classic FIR filter');
 attenuation=ridgepack_db(0.001);
 ripple=ridgepack_db(0.98);
 if nargin<6 | isempty(bandwidth) ; bandwidth=0.5; end
 firc=true;

elseif strcmp(kind,'FIRH') 

 disp('Hibler FIR filter');
 attenuation=ridgepack_db(0.001);
 ripple=ridgepack_db(0.99);
 if nargin<6 | isempty(bandwidth) ; bandwidth=0.5; end
 firh=true;
 if ~strcmp(type,'low')
  error('For the Hibler cosine filter, only a low pass is possible in this program')
 end

elseif strcmp(kind,'IIR') 

 disp('IIR Butterworth filter');
 attenuation=ridgepack_db(0.05);
 ripple=ridgepack_db(0.95);
 if nargin<6 | isempty(bandwidth) ; bandwidth=1.0; end
 iirb=true;

else

 error('Filter type is incorrect')

end

if nargin<7; smooth=[] ; end

if any(Wn>fs/2)
	error(['nth elements of Wn greater than the Nyquist frequency: ',...
	       num2str(find(Wn>fs/2))])
elseif bandwidth>=fs/2
	disp('WARNING: Bandwidth greater than Nyquist frequency')
	bandwidth=4*fs/10;
	disp(['Reseting bandwidth to ',num2str(bandwidth),' cycles/day'])
end

% Determine the search range, or else set n to the specified even order.
% Set the minimum attenuation and ripple criteria if n is directory specified
if emptyn
 if firc | firh
	nmin=20; nmax=2000;
 elseif iirb
	nmin=2; nmax=50;
 else
	error('Incorrect filter settings')
 end
 disp(['Bandwidth tolerance is:   ',num2str(bandwidth)]);
 disp(['Attenuation tolerance is: ',num2str(ridgepack_invdb(attenuation))]);
 disp(['Ripple tolerance is:      ',num2str(ridgepack_invdb(ripple))]);
 disp('-------------------------------');
else
 n=2*ceil(n/2); nmin=n; nmax=n;
end

% find the minimum even filter order for the given specifications
for n=nmin:2:nmax;

  % diagnostic
  if floor(n/100)==n/100
	disp(['Still searching at filter order ',num2str(n)]);
  end

  % Get difference equation coefficients
  if firc
  	a=1; b=fir1(n,Wn*2/fs,type,hamming(n+1),'scale');
  elseif firh
  	a=1; b=ridgepack_fir1hibler(n,Wn*2/fs);
  elseif iirb
  	[b,a] = butter(n,Wn*2/fs,type);
  else
	error('Incorrect filter settings')
  end

  % get frequency response
  [H,W]=freqz(b,a,n); W=(fs/2)*(W/pi);

  % determine if this frequency response passes criteria and pass by some checks
  % if the filter order has been specified.
  if ~emptyn
   pass=true;
   description{1}='Filter';
  else
   [pass,description]=ridgepack_fcriteria(H,W,Wn,bandwidth,attenuation,ripple,ts.timeinfo.unit);
  end

  % if nmax is reached, exit with an error
  if n==nmax & nmax==nmin & ~pass; 
   disp(description);
   error(['Filter did not meet specifications up to order ',num2str(n)]); 
  elseif n==nmax & ~pass; 
   disp(description);
   error(['Filter did not meet specifications up to order ',num2str(n)]); 
  end
 
  % if criteria are met, exit from the loop
  if pass; break; end

end

% display
disp(['Filter order is ',num2str(n)]);

% add filter type to description
if firc
  descrip{1}=['FIR ',char(description{1})];
  description=descrip;
elseif firh
  descrip{1}=['FIR Hibler ',char(description{1})];
  description=descrip;
elseif iirb
  descrip{1}=['IIR ',char(description{1})];
  description=descrip;
end


% zero-phase digital filtering with minimized startup transients
if nargout>0

	oldsize=size(ts.Data);
	newsize=[oldsize(1), prod(oldsize(2:end))];
	ts.Data=reshape(ts.Data,newsize);

        wavelength=1/Wn(1);
        samplength=1/fs;
        nlength=ceil(wavelength/samplength);

	disp(['Filtering sequence of ',num2str(newsize(2)),' timeseries...']);
        	
	% Parallelize this loop if possible
	for i=1:newsize(2)


         if ~isempty(smooth) & (strcmp(type,'low') | strcmp(type,'pass'))
          ts.Data(1:nlength,i)=smooth+(ts.Data(nlength+1,i)-smooth)*...
                               sin(2*pi*[1:nlength]*samplength/(4*wavelength))
         end
  
         if length(ts.Data(:,i))<3*(length(b)-1)
	  error(['Filter order ',num2str(length(b)-1),...
                 ' too large for data length ',num2str(length(ts.Data(:,i)))])
         end

	 ts.Data(:,i)=filtfilt(b,a,ts.Data(:,i));

	end

	disp('done.')

	ts.Data=reshape(ts.Data,oldsize);

else

	% plot filter characteristics on a graph
	clf;
	ncfreqzplot(H,W,b,a,description,ts.timeinfo.unit);
	clear ts;

end


disp('-------------------------------');
disp(['Maximum frequency response: ',num2str(max(abs(H(:))))])
disp(['Minimum frequency response: ',num2str(min(abs(H(:))))])
disp('-------------------------------');

end % function
