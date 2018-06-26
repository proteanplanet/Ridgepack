function [crosscov]=ridgepack_cov(ts1,ts2,maxlags)

% ridgepack_cov - User defined biased cross covariance estimator
%
% function [crosscov]=ridgepack_cov(ts1,ts2,maxlags)
%
% This function calculates the cross covariance of two timeseries, 
% ts1 and ts2, which are matlab timeseries objects.  The cross
% covariance is approximated for a maximum lag of maxlags, and
% is a biased estimator, meaning that the covariance for each 
% sample is divided by the length, N, of the timeseries, rather
% than N-m where m is the lag.  For more information, see 
% Priestley (1981), "Spectral Analysis and Time Series".
%
% INPUT:
%
% ts1 - Matlab time series object
% ts2 - Matlab time series object of same length and time samples as ts1
% maxlags - maximum length of the cross covariance series
%
%
% OUTPUT:
%
% crosscov - cross covariance series for ts1 and ts2
%
% Note that the timeseries can be either real or complex.
%
% Ridgepack Version 1.0
% Andrew Roberts, Naval Postgraduate School, March 2018 (afrobert@nps.edu)
% Naval Postgraduate School, Department of Oceanography
%

if (ts1.Time~=ts2.Time)
	error('Time series are not synchronized')
end

N=length(ts1.Data);
m1=mean(ts1.Data);
m2=mean(ts2.Data);

m=zeros([(2*N)-1 1]);
crosscov=zeros([(2*N)-1 1]);

for i=1:N
for j=1:N

 k=i-j+N;

 if(abs(i-j)<=maxlags) | maxlags==0

  m(k)=m(k)+1;

  crosscov(k)=crosscov(k)+(ts1.Data(i)-m1)*conj(ts2.Data(j)-m2);

 end

end
end

crosscov=crosscov(m>0)./N;


