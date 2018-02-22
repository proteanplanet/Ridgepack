function idx=ridgepack_nearestval(x,y)

% ridgepack_nearestval - provides an array's index of nearest value to that specified
%
% function idx=nearestval(x,y)
%
% This provides the index of the number nearest to that specified
% for matrix x and value for which you are searching=y.
%
% This is simply the matlab line:
% idx = find(abs(x(:)-y)==min(abs(x(:)-y)));
%
% Andrew Roberts, Naval Postgraduate School, March 2018  (afrobert@nps.edu)
%

global debug;
if debug; disp(['Entering ',mfilename,'...']); end

idx = find(abs(x(:)-y)==min(abs(x(:)-y)));

if debug; disp(['...Leaving ',mfilename]); end

