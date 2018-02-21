function idx=ncnearestval(x,y)

% NCNEARESTVAL - provides an array's index of nearest value to that specified
%
% function idx=nearestval(x,y)
%
% This provides the index of the number nearest to that specified
% for matrix x and value for which you are searching=y.
%
% This is simply the matlab line:
% idx = find(abs(x(:)-y)==min(abs(x(:)-y)));
%
% Written by Andrew Roberts 2012
% Department of Oceanography, Naval Postgraduate School
%
% $Id$  

global debug;
if debug; disp(['Entering ',mfilename,'...']); end

idx = find(abs(x(:)-y)==min(abs(x(:)-y)));

if debug; disp(['...Leaving ',mfilename]); end

