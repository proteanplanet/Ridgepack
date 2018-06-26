function [b]=fir1hibler(n,Wp)

% FIR1HIBLER - Hibler 1972 FIR time series filter
%
% function [b]=fir1hibler(n,Wp)
%
% Low pass FIR filter following Hibler, 1972.
% 
% This function behaves in a similar manner to the generic fir1 
% matlab routine.  However the filter is slightly different, 
% constructed of cosine filter weights and filtered as a convolution.
% This provides the filter weights for the b component of the filter
% weights for an FIR filter (a=1).  For more information see fir1. 
%
% INPUT:
%
% n	- Order of the filter.  This must be an even number.
% Wp    - Maximum passband frequency of the filter between 0 and 1, 
%         where 1 corresponds to the Nyquist frequency. For this
%         particular filter, this corresponds to N1/(2*N) cycles
%         per data interval, with the transition band extending 
%         to (N1+3)/(2*N). Ripple errors outside the transition
%         band are <0.9%. Note that N1 is determined by rounding
%         to the nearest integer value using Wp.
%
%
% OUTPUT:
%
% b	- A row vector b containing the n+1 coefficient of an 
%         order n lowpass FIR filter.
%
% For more information on this filter design, see the Appendix of:
% Hibler, W. D., III, 1972: Removal of Aircraft motion from LASER 
% profilometers. J. Geophys. Res., 77, 36.
% 
% Ridgepack Version 1.0
% Andrew Roberts, Naval Postgraduate School, March 2018 (afrobert@nps.edu)
% Naval Postgraduate School, Department of Oceanography
%

% Add 1 to the order for n+1 filter coefficients
n=n+1;

% Make sure the n+1 is an odd number
if fix((n-1)/2) ~= (n-1)/2
 error('The order of the filter must be even')
end

% set up N for the filter code based on the order
N=(n-1)/2;

% set N1 accordingly for the edge of the passband
N1=round(2*N*Wp);

disp(['Passband ends at ',num2str(N1/(2*N)),', stopband starts at ',num2str((N1+3)/(2*N))]);

% Allocate C(n) and H(n)
C=ones((2*N)+1, 1);
H=ones(N+1, 1);

% Set coefficients
for i=0:N
 if i<=N1
  H(i+1)=1;
 elseif i==N1+1
  H(i+1)=0.77;
 elseif i==N1+2
  H(i+1)=0.23;
 else
  H(i+1)=0.;
 end
end

for n=-N:N

 C(n+N+1)=0;

 for i=1:N-1

  C(n+N+1)=C(n+N+1)+H(i+1)*cos(pi*n*i/N)/N;

 end

 C(n+N+1)=C(n+N+1)+H(1)/(2*N)+cos(pi*n)*H(N+1)/(2*N);

end

C(1)=C(1)/2;

C(end)=C(end)/2;

if sum(C(1:end)) ~= 1
 disp(['Filter weights do not sum to one: sum=',num2str(sum(C(1:end)))]);
end

b=C;

end % function


