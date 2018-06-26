function [hd]=ridgepack_ideallpf(wc,M)

% ridgepack_ideallpf - Ideal low pass filter computation
%
% function hd=ridgepack_ideallpf(wc,M)
%
% This is an ideal low pass filter computation.
%
% INPUT:
%
% wc = cutoff frequency in radians
% M = length of the ideal filter
%
%
% OUTPUT:
%
% hd = ideal impulse response between 0 to M-1
%
% Ridgepack Version 1.0
% Andrew Roberts, Naval Postgraduate School, March 2018 (afrobert@nps.edu)
% Naval Postgraduate School, Department of Oceanography
%

alpha = (M-1)/2;
n=[0:1:(M-1)];
m=n-alpha+eps;
hd=sin(wc*m)./(pi*m);

end

