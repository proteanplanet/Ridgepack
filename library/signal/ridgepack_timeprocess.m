function [nc1,nc2]=ridgepack_timeprocess(nc1,name1,nc2,name2,interp,resample,Wn,firl,type,n,bandwidth,smooth)

% ridgepack_timeprocess - Interpolates and filters a time series in an nc structure
%
% function [nc1,nc2]=ridgepack_timeprocess(nc1,name1,nc2,name2,interp,resample,Wn,firl,type,n,bandwidth,smooth)
%
% This function combines the ridgepack_2timeseries/ridgepack_2vectseries functions with ridgepack_filter1
% to process a time series in an nc structure and provide output in an nc structure
% that is resampled and filtered.  Inputs follow ridgepack_2vectseries and ridgepack_filter1.
% If no filter cutoff frequencies are specified with Wn, then no filtering is performed
% on the timeseries.  In other words, all arguments after name2 are optional, and
% if all arguments after and including Wn are omitted, no filtering is performed,
% only resampling.  See manual pages for the functions ridgepack_2timeseries, ridgepack_2vectseries
% and ridgepack_filter1 for further information on the inputs.
%
% Ridgepack Version 1.0
% Andrew Roberts, Naval Postgraduate School, March 2018 (afrobert@nps.edu)
% Naval Postgraduate School, Department of Oceanography

global debug;

if nargin<10; n=[]; end
if nargin<11; bandwidth=[]; end

if debug; disp('Resampling timeseries...'); end

if nargin<6; error('Must specify interpolation type and reampling'); end

if isempty(nc2) | strcmp(nc2,'')
 ts=ridgepack_2timeseries(nc1,name1,interp,resample);
else
 if isempty(setdiff(inputname(1),inputname(3)))
  ts=ridgepack_2vectseries(nc1,name1,nc1,name2,interp,resample);
 else
  ts=ridgepack_2vectseries(nc1,name1,nc2,name2,interp,resample);
 end
end

if debug; disp('Filtering timeseries...'); end

if nargin<7; 
 disp('No filtering is being performed on this time series')
elseif nargin<8; 
 [ts,b,a,description]=ridgepack_filter1(ts,Wn);
elseif nargin<9; 
 [ts,b,a,description]=ridgepack_filter1(ts,Wn,firl);
elseif nargin<10; 
 [ts,b,a,description]=ridgepack_filter1(ts,Wn,firl,type);
elseif nargin<11; 
 [ts,b,a,description]=ridgepack_filter1(ts,Wn,firl,type,n);
elseif nargin<12; 
 [ts,b,a,description]=ridgepack_filter1(ts,Wn,firl,type,n,bandwidth);
else
 [ts,b,a,description]=ridgepack_filter1(ts,Wn,firl,type,n,bandwidth,smooth);
end

if debug; disp('Adding metadata...'); end

nc1.time.data=ts.Time;
nc1.(name1).data=real(ts.Data);
if nargin<7; 
 nc1.(name1).long_name=['Resampled ',nc1.(name1).long_name];
else
 nc1.(name1).long_name=['Filtered ',nc1.(name1).long_name];
end

if isfield(nc1.attributes,'comment') & nargin>6
   nc1.attributes.comment=[nc1.attributes.comment,' Temporally filtered using ',...
   nccellcat(description)];
elseif nargin>6
   nc1.attributes.comment=['Temporally filtered using ',nccellcat(description)];
end

if ~isempty(nc2);

 if isempty(setdiff(inputname(1),inputname(3)))

  nc1.(name2).data=imag(ts.Data);
  if nargin<7; 
   nc1.(name2).long_name=['Resampled ',nc1.(name2).long_name];
  else
   nc1.(name2).long_name=['Filtered ',nc1.(name2).long_name];
  end

  if isfield(nc1.attributes,'comment') & nargin>6
   nc1.attributes.comment=[nc1.attributes.comment,' Time filtered using ',...
   nccellcat(description)];
  elseif nargin>6
   nc1.attributes.comment=['Time filtered using ',nccellcat(description)];
  end

 else

  nc2.time.data=ts.Time;
  nc2.(name2).data=imag(ts.Data);
  if nargin<7; 
   nc2.(name2).long_name=['Resampled ',nc2.(name2).long_name];
  else
   nc2.(name2).long_name=['Filtered ',nc2.(name2).long_name];
  end

  if isfield(nc2.attributes,'comment') & nargin>6
   nc2.attributes.comment=[nc2.attributes.comment,' Time filtered using ',...
   nccellcat(description)];
  elseif nargin>6
   nc2.attributes.comment=['Time filtered using ',nccellcat(description)];
  end

 end

end


