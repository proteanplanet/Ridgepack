function [nc,ts]=ridgepack_waveform(frequency,phase,tlength,noise,amp,interval);

% ridgepack_waveform - Generates a 1D complex waveform time series
%
% function [nc,ts]=ridgepack_waveform(frequency,phase,tlength,noise,amp,interval);
%
% This function generates a 1D sinusoid complex waveform as a time series.
% The waveform contains both clockwise and anticlockwise components.  To
% generate clockwise, enter negative frequencies, to generate anticlockwise
% timeseries, enter positive frequencies.
% 
% INPUT:
%
% frequency - Wavelength in terms of cycles per day. This may be a vector
%              if two or more frequencies are required. 
% phase     - Phase of the wave form as a fraction from 0 to 1 of the period
%             The phase must correspond be a vector corresponding to each 
%             frequency provided.
% tlength   - Number of days for which the timeseries is to last
% noise     - If set to 'true' third arguement is provided, white noise
%             random numbers) is added to the timeseries from a normal 
%             distribution.
% amp       - Optional fourth argument giving the amplitude for each
%             specified frequency (1 is the default)
% interval  - Sampling interval in minutes (default is 10 minutes).
%
%
% OUTPUT:
%
% nc - nc structure containing the waveform timeseries tlength days long	
% ts - matlab timeseries structure
%
% The wave forms have a default sampling period of 10 minutes.
%
% Ridgepack Version 1.0
% Andrew Roberts, Naval Postgraduate School, March 2018 (afrobert@nps.edu)
% Naval Postgraduate School, Department of Oceanography
%

if nargin<3; error('Insufficient number of arguments'); end
if nargin<5; 
 amp=ones(size(frequency))  
elseif length(amp)~=length(frequency)
 error('amplitude must have the same number of elements as frequency')
end

if nargin<6
 interval=10;
elseif interval>1440 | interval<1
 error('Interval must be less than 1 day and greater than 1 minutes')
end

% write title
nc.attributes.title=['Time series waveform'];

% generate time series with 10 minute intervals
nc.time.data=[0:interval/1440:tlength];
nc.time.long_name='time';
nc.time.units='days since 0000-00-00 00:00:00';
nc.time.dimension={'time'};
nc.time.type='NC_DOUBLE';

% generate waveform, beinc careful to not simply add complex conjugates together to negate
% the clockwise rotation sense.
wave=zeros(length(nc.time.data),length(frequency));
for i=1:length(frequency)
 if frequency(i)>0 % anticlockwise
  wave(:,i)=amp(i)*cos(2*pi*frequency(i)*nc.time.data+phase(i)*2*pi)+amp(i)*sqrt(-1)*sin(2*pi*frequency(i)*nc.time.data+phase(i)*2*pi);
 else % clockwise
  wave(:,i)=amp(i)*cos(2*pi*abs(frequency(i))*nc.time.data+phase(i)*2*pi)-amp(i)*sqrt(-1)*sin(2*pi*abs(frequency(i))*nc.time.data+phase(i)*2*pi);
 end
end
nc.wave.data=sum(wave,2);
nc.wave.long_name=['Complex sinusoid waveform timeseries with frequencies: ',num2str(frequency),' days'];
nc.wave.dimension={'time'};
nc.wave.type='NC_DOUBLE';

% add white noise to waveform if required
if noise; nc.wave.data=nc.wave.data+random('Normal',0,1,length(nc.wave.data),1); end

% arrange structure
nc=ncstruct(nc);

% write the timeseries to a timeseries structure
ts=timeseries(squeeze(nc.wave.data),nc.time.data,'Name',nc.wave.long_name);

% write data information
ts.DataInfo.UserData=nc.attributes.title;

% write time information
ts.TimeInfo.Units='days';
ts.TimeInfo.Format='dd-mmm-yyyy HH:MM:SS';

% set time increment
ts.TimeInfo.Increment=datenum(ts.Time(2)-ts.Time(1));

end

