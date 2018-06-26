function [ts]=ridgepack_2vectseries(nc1,name1,nc2,name2,interp,resample)

% ridgepack_2vectseries - Creates a complex timeseries object for a 2D vector 
%
% function [ts]=ridgepack_2vectseries(nc1,name1,nc2,name2,interp,resample)
%
% This function creates a timeseries object, ts, given:
%
% INPUT:
%
% nc1      - nc structure in which the first time series is located
% 	     (see ncstruct for more information)
% name1    - string providing the variable to be converted into 
%            a timeseries for the u-component of a vector.
% nc2      - nc structure in which second time series is located
%            if the time series is a 2D vector such as velocity.
% name2    - string providing the second variable to be used in 
%            converting the timeseries into 2D vector timeseries.
%            This variable should have the same time vector and
%            be orthogonal to nc1.name1, being the v-component
%            of the vector.
% interp   - method of interpolation (optional).  This many be
%            set to 'linear', 'spline', 'cubic' or other 
%            interpolation methods listed under the Matlab 
%            function "setinterpmethod".  The default is 'linear'.
% resample - [year month day hour minute second.fraction] vector
%            specifying the new time increment for resampling
%            the time series.
%
%
% OUTPUT:
%
% ts    - Matlab timeseries object
%
% A 2D velocity timeseries object may be created by
% including the names of two fields, name1 and name2. For
% example, name1='u' and name2='v'. The timeseries is stored
% as a complex number ts.Data=u+sqrt(-1)*v. A fourier transform
% or Power Spectral Density taken of the timeseries will then 
% automatically be rotary, or two sided, without any additional
% information needing to be provided by the user. 
%
% This function synchronizes the two timeseries before adding 
% them together.
%
% Ridgepack Version 1.0
% Andrew Roberts, Naval Postgraduate School, March 2018 (afrobert@nps.edu)
% Naval Postgraduate School, Department of Oceanography
%

% check correct number of arguments as input
if nargin<2; error('Insufficient input arguments'); end
if nargin==2; nc2=''; name2=''; end
if nargin==3; error('Must specify name of second input field'); end
if nargin < 5; interp='linear'; end
if nargin<6; resample=[]; end

% check that nc1 and nc2 are structures
if ~isstruct(nc1); error(['Input nc1 is not a structure.']); end
if ~strcmp(nc2,'') & not(isstruct(nc2)); error(['Input nc2 is not a structure.']); end

% if nc1 and nc2 are identical, then there is no need for synchronization 
if isempty(setdiff({inputname(1)},{inputname(3)}))

 disp('Processing dual timeseries as a single complex number')

 nc1.(name1).data=nc1.(name1).data+sqrt(-1)*nc2.(name2).data;

 ts=ridgepack_2timeseries(nc1,name1,interp,resample);

else

 disp('Synchronizing dual timeseries')

 % create independent time series
 ts1=ridgepack_2timeseries(nc1,name1,interp,resample);
 ts2=ridgepack_2timeseries(nc2,name2,interp,resample);

 % make sure they are synchronized
 %if isnan(ts1.TimeInfo.Increment)
 % error('ts1 does not have a time increment set')
 %else
 % deltat=ts1.TimeInfo.Increment;
  deltat=mean(ts1.Time(2:end)-ts1.Time(1:end-1));
  [ts1 ts2]=synchronize(ts1,ts2,'uniform','interval',deltat);
 %end

 % Now convert to a single timeseries
 ts=ts1;
 ts.Data=ts1.Data+sqrt(-1)*ts2.Data;
 ts.TimeInfo.Increment=deltat;

end

