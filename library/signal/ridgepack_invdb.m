function [power]=ridgepack_invdb(decibels)

% ridgepack_invdb - Convert decibels to power
%
% function [power]=ridgepack_invdb(decibels)
%
% INPUT:
% 
% decibels - value converted to decibels using 10*log10(power)
%
%
% OUTPUT:
%
% Power - linear value to be converted from decibels that is not in units of volts.
%
% Ridgepack Version 1.0
% Andrew Roberts, Naval Postgraduate School, March 2018 (afrobert@nps.edu)
% Naval Postgraduate School, Department of Oceanography
%

power=10.^(decibels/10);

end % function

