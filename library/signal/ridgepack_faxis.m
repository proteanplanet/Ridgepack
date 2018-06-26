function [faxis,fs,df]=ridgepack_faxis(ts)

% ridgepack_faxis - Returns a shifted frequency axis, sampling frequency and frequency resolution
%
% function [faxis,fs,df]=ridgepack_faxis(ts)
%
% Returns a frequency axis, Nyquist frequency and frequency resolution 
% given a time series object.  The frequency axis is shifted from the range
% [0,fs) to (-fs,fs] if length(ts) is even, or else (-fs,fs) if length(ts)
% is odd.  All frequency axes include the frequency of zero.
%
% INPUT:
%
% ts - timeseries object generated with ridgepack_2timeseries
%
%
% OUTPUT:
%
% faxis - frequency axis of fourier transform of ts shifted with fftshift
% fs    - 2*Nyquist frequency = sampling frequency
% df    - frequency resolution 
%
% Ridgepack Version 1.0
% Andrew Roberts, Naval Postgraduate School, March 2018 (afrobert@nps.edu)
% Naval Postgraduate School, Department of Oceanography
%

% get length of fft
N=length(ts.data);

% get the sampling period
if isnan(ts.TimeInfo.Increment) & isempty(ts.TimeInfo.UserData);
 error('Samples in time series are not uniformly sampled')
elseif isnan(ts.TimeInfo.Increment) & ischar(ts.TimeInfo.UserData) & ...
       length(ts.TimeInfo.UserData)>3 & strcmpi(ts.Timeinfo.UserData(1:3),'dt=')
 try
  tsamp=str2num(ts.TimeInfo.UserData(4:end));
 catch
  error('Unable to convert ts.Timeinfo.UserData to a time interval')
 end
elseif isempty(ts.TimeInfo.UserData);
 tsamp=ts.TimeInfo.Increment;
elseif strcmpi(ts.TimeInfo.UserData(1:3),'dt=')
 tsamp=str2num(ts.TimeInfo.UserData(4:end));
else
 error('Samples in time series are not uniformly sampled')
end

% get sampling frequency and frequency resolution
fs=1/tsamp; df=fs/N;

if floor(N/2)==N/2 % even case
	faxis=[1-(N/2):1:(N/2)]*df;
else % odd case
	faxis=[(1-N)/2:1:(N-1)/2]*df;
end

faxis=faxis';

% check the output for stupid coding errors 
if isinf(fs) | isnan(fs)
	error('fn is either infinite or not a number')
elseif length(faxis)~=N
	disp(num2str(length(faxis)))
	disp(num2str(N))
	error('frequency axis is not correct length')
end

