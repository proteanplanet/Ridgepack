function [ts,deltat]=ridgepack_uniform_sample_check(ts)

% ridgepack_uniform_sample_check - Checks for a uniform sampling period in timeseries objects
% 
% function [ts,deltat]=ridgepack_uniform_sample_check(ts)
% 
% This function checks for a uniform sampling period in a Matlab timeseries 
% object and if none is found, checks to see if there is only a minor inaccuracy
% in the sampling period.  If this is the case, the sampling times are rounded without
% changing the values of the timeseries, else an error is registered suggesting 
% the the timeseries is resampled.
%
% INPUTS:
% ts - matlab time series object
%
%
% OUTPUT:
%
% deltat - timestep of timeseries
%
% Ridgepack Version 1.0
% Andrew Roberts, Naval Postgraduate School, March 2018 (afrobert@nps.edu)
% Naval Postgraduate School, Department of Oceanography
%


if isnan(ts.TimeInfo.Increment) 
 deltat=mean(ts.Time(2:end)-ts.Time(1:end-1));

 % check for accuracy in timestep to six significant digits
 md=6-floor(log10(max(diff(ts.Time))));
 err=max(diff(round(diff(ts.Time)*10^md)));
 err=max(diff(fix(diff(ts.Time)*10^md)));

 % if accurate to six decimal places, round the sampling to the nearest 1000th of second
 if err==0;
  disp(['Samples not perfectly separated in time, mean interval is: ', num2str(mean(deltat))]);
  dt=max(round(diff(ts.Time)*10^md))*(10^-md);
  dt=max(fix(diff(ts.Time)*10^md))*(10^-md);
  resample=datevec(dt);
  if resample(6)<10^-3 & any(resample(1:5)>0)
   resample(6)=0.0;
  elseif (60-resample(6))<10^-3 & any(resample(1:5)>0)
   resample(5)=resample(5)+1;
   resample(6)=0.0;
  end
  disp(['Linearly resampling to a timestep of exactly ',datestr(datenum(resample),'dd+HH:MM:SS')])
 else
  disp(['Samples not uniformly distributed in time, mean interval is: ', num2str(mean(deltat))]);
  error('Please use the resample facility in this function');
 end
else
 deltat=ts.TimeInfo.Increment;
end


