function [ts]=ridgepack_2timeseries(nc1,name1,interp,resample)

% ridgepack_2timeseries - Creates a timeseries object from an nc1 stucture variable
%
% function [ts]=ridgepack_2timeseries(nc1,name1,interp,resample)
%
% This function creates a timeseries object, ts, given:
%
% INPUT:
%
% nc1      - nc1 structure in which first time series is located
% 	     (see ncstruct for more information)
% name1    - string providing the variable to be converted into 
%            a timeseries.
% interp   - method of interpolation (optional).  This many be
%            set to 'linear', 'spline', 'cubic', 'pchip' or other 
%            interpolation methods listed under the Matlab 
%            function "setinterpmethod".  The default is 'linear'.
%            The interp method can be entered as 'constant' if
%            the first time step can be assumed to be the constant
%            timestep throughout the timeseries.  If this is the case
%            no resample value needs to be specified.
% resample - [year month day hour minute second.fraction] vector
%            specifying the new time increment for re-sampling
%            the time series. Hourly is assumed if nothing is inserted.
%
%
% OUTPUT:
%
% ts    - Matlab timeseries
%
% Ridgepack Version 1.0
% Andrew Roberts, Naval Postgraduate School, March 2018 (afrobert@nps.edu)


% check correct number of arguments as input
if nargin<2; error('Insufficient input arguments'); end
if nargin<3; interp='linear'; end
if nargin<4; 
 resample=[0 0 0 1 0 0.0]; 
elseif ~isnumeric(resample)
 error('resample must be numeric with format [year month day hour minute second])')
elseif length(resample)~=6
 error('resample must have format [year month day hour minute second])')
end

% check that nc1 is a structure
if not(isstruct(nc1)); error([inputname(1),' is not a structure']); end

% Position time as first dimension of the variable, making 
% sure that time is indeed a dimension of the variable.
if isempty(intersect({'time'},nc1.(name1).dimension))
 error([name1,' does not have time as a dimension']);
else
 otherdims=ridgepack_setdiff(nc1.(name1).dimension,{'time'});
 if isempty(otherdims)
  nc1.(name1).data=nc1.(name1).data(:);  
 else
  nc1=ridgepack_shuffle(nc1,otherdims);
 end
end

% Get data and standard units
[z,nc1]=ridgepack_standardunits(nc1,name1);
clear z;

% Remove general NaNs from data if dealing with 
% 2D array only (including time) and remove samples
% where there are large gaps
if any(isnan(nc1.(name1).data(:))) & ndims(nc1.(name1).data)==2
 disp('Removing NaN tainted records from the 2D array');
 j=0;
 for i=1:length(nc1.time.data)
	if max(isnan(nc1.(name1).data(i,:)))<1
		j=j+1;
		datx(j,:)=nc1.(name1).data(i,:);
		time(j)=nc1.time.data(i);
	end
 end

 if isvarname('time')

  tdiff=NaN;
  time=nc1.time.data;
  datx=NaN*ones(size(nc1.(name1).data));

 else

  % add NaNs where there are large gaps in the data:
  % find large gaps in timeseries and remove data outside three
  % standard deviations great than the normal time difference
  tdiff=time(2:length(time))-time(1:length(time)-1);
  stdtdiff=median(tdiff);
  disp(['Adding NaNs at time gaps greater than ',num2str(10*stdtdiff)]);
  for i=1:length(time)-1
	if tdiff(i)>10*stdtdiff
		datx(i,:)=NaN;
		datx(i+1,:)=NaN;
		disp(['...adding NaN at ',num2str(i)])
	end
  end

 end

else
 datx=nc1.(name1).data;
 time=nc1.time.data;
end

% write data, time and name
ts=timeseries(datx,time,'Name',nc1.(name1).long_name);

% treat NaN as missing value
ts.TreatNaNasMissing=true;

% write data information
if isfield(nc1.(name1),'units')
 ts.DataInfo.Unit=nc1.(name1).units;
end
ts.DataInfo.UserData=nc1.attributes.title;

% set interpolation method
myFuncHandle=@(new_Time,Time,Data) interp1(Time,Data,new_Time,interp,NaN);
myInterpObj=tsdata.interpolation(myFuncHandle);
ts=setinterpmethod(ts, myInterpObj);

% write time information
ts.TimeInfo.Units='days';
ts.TimeInfo.Format='dd-mmm-yyyy HH:MM:SS';

% set constant time increment if requested
if strcmp('constant',interp)
 ts.TimeInfo.Increment=datenum(ts.Time(2)-ts.Time(1));
end

% check that all time intervals are identical, and resample if not
if isempty(resample)
 [ts,deltat]=ridgepack_uniform_sample_check(ts);
end

% resample if required, using hours as the base unit
if ~isempty(resample)
 disp('Resampling the timeseries')
 deltat=datenum(resample);
 newtime=ts.Time(1):deltat:ts.Time(end);
 ts=ridgepack_resample(ts,newtime,interp);
 ts.Timeinfo.UserData=['dt=',num2str(deltat)];
end


