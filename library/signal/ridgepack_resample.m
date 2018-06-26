function [tsnew]=ridgepack_resample(ts,newtime,method)

% ridgepack_resample - Resample time series using 1D spline interpolation (no extrapolation)
%
% function [tsnew]=ridgepack_resample(ts,newtime,method)
%
% This function directly replaces the interpolation object in a time series ts
% with another interpolation method without extrapolation.  Where extrapolation 
% would be required a NaN is placed in the output time series which is resampled to the 
% supplied time vector. The output is tsnew, a Matlab time series object.
%
% INPUTS:
% ts      - input time series
% newtime - input new time vector to which the ts is to be resampled
% method  - options are 'linear' or 'spline'.  If omitted, the default
%           is 'linear'.  This option can also be set to 'extrap' for
%           linear interpolation with extrapolation.
%
% OUTPUTS:
% tsnew   - output time series that is re-sampled
%
% Ridgepack Version 1.0
% Andrew Roberts, Naval Postgraduate School, March 2018 (afrobert@nps.edu)
% Naval Postgraduate School, Department of Oceanography
%

% check for arguments
if nargin<1
 error('No time series object supplied')
end

% set defaults
if nargin<3
 method='linear';
end

tsnew=ts;

if strcmp(method,'extrap')
 % extrapolation
 myFuncHandle = @(new_Time,Time,Data) interp1(Time,Data,new_Time,'linear','extrap');
elseif strcmp(method,'pchip')
 % spline with no extrapolation - fill NaNs
 myFuncHandle = @(new_Time,Time,Data) interp1(Time,Data,new_Time,'pchip',NaN);
elseif strcmp(method,'spline')
 % spline with no extrapolation - fill NaNs
 myFuncHandle = @(new_Time,Time,Data) interp1(Time,Data,new_Time,'spline',NaN);
elseif strcmp(method,'linear')
 % linear no extrapolation - fill NaNs
 myFuncHandle = @(new_Time,Time,Data) interp1(Time,Data,new_Time,'linear',NaN);
else
 error('method chosen not usable')
end

myInterpObj = tsdata.interpolation(myFuncHandle);
tsnew = setinterpmethod(tsnew, myInterpObj);
tsnew = resample(tsnew,newtime);

tsnew.TimeInfo

disp(['Resampled using ',method,' method']);

tsnew.Name=['Resampled from ',ts.Name];

tsnew.DataInfo.UserData=['Sourced from ',tsnew.DataInfo.UserData];

