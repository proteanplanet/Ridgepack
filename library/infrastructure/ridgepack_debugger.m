function ridgepack_debugger(signal)

% ridgepack_debugger - Turn degging off and on in icepack
%
% function ridgepack_debugger(signal)
%
% This function turns on or off debugging for the entire MATLAB
% icepack package.
%
% INPUT:
%
% signal - Character string, either 'on' or 'off'.  Any other 
%          string is rejected.
%
% Ridgepack Version 1.0
% Andrew Roberts, Naval Postgraduate School, March 2018 (afrobert@nps.edu)

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
 
