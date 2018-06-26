function [decibels]=ridgepack_db(power)

% ridgepack_db - Convert power to decibels
%
% function [decibels]=ridgepack_db(power)
%
% INPUT:
% 
% Power - linear value to be converted to decibels that is not in units of volts.
%
%
% OUTPUT:
%
% decibels - value converted to decibels using 10*log10(power)
%
% Ridgepack Version 1.0
% Andrew Roberts, Naval Postgraduate School, March 2018 (afrobert@nps.edu)
% Naval Postgraduate School, Department of Oceanography
%

decibels=10*log10(power);

