function [xk]=ridgepack_dft(xn)

% ridgepack_dft - user defined Discrete Fourier Transform
%
% function [xk]=ridgepack_dft(xn)
%
% This is a user defined Discrete Fourier Transform.
% For practice purposes.  Use fft when possible.
%
% INPUT:
%
% xn - N-point finite-duration sequence
%
%
% OUTPUT:
%
% Xk - DFT coefficient array 0<=k<=N-1
%
% Ridgepack Version 1.0
% Andrew Roberts, Naval Postgraduate School, March 2018 (afrobert@nps.edu)
% Naval Postgraduate School, Department of Oceanography
%


N=length(xn);
n=[0:1:N-1]';
k=[0:1:N-1];
WN=exp(-j*2*pi/N);
nk=n*k;
WNnk=WN.^nk;
xk=WNnk*xn;


