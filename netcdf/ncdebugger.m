function ncdebugger(signal)

% NCDEBUG - Turn degging off and on in icepack
%
% function ncdebug(signal)
%
% This function turns on or off debugging for the entire MATLAB
% icepack package.
%
% Input:
% signal - Character string, either 'on' or 'off'.  Any other 
%          string is rejected.
%
% Written by Andrew Roberts, October 2012
% Naval Postgraduate School, Department of Oceanography
%
% $Id$

global debug;

if nargin<1
 error('ncdebug missing input signal')
elseif ~ischar(signal)
 error('Input signal must be a character string ''on'' of ''off''')
elseif strcmp(signal,'on')
 debug=1;
elseif strcmp(signal,'off')
 debug=0;
else
 error('Input signal must be either ''on'' of ''off''')
end

disp(['Debugger is ',signal])
 
