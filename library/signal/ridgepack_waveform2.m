function nc=ridgepack_waveform2(nc,wavelength);

% ridgepack_waveform2 - Generates a 2D waveform on a given grid
%
% function nc=ridgepack_waveform2(nc,wavelength);
%
% This function generates a 2D sinusoid waveform on a given grid 
% provided in the input nc structure.
% 
% INPUT:
%
% nc         - nc structure from icepack with an cartesian grid
% wavelength - wavelength in terms of grid points of waveform
%
% Ridgepack Version 1.0
% Andrew Roberts, Naval Postgraduate School, March 2018 (afrobert@nps.edu)
% Naval Postgraduate School, Department of Oceanography
%

if ~isstruct(nc);
	error('nc is not a structure');
end

nc=ncshuffle(nc,{'x','y'});

[zx,zy]=meshgrid(0.5*sin(nc.x.data*2*pi/wavelength),0.5*sin(nc.y.data*2*pi/wavelength));

nc.zwave.data=zx+zy;
nc.zwave.long_name=['Sinusoid waveform with wavelength ',num2str(wavelength)];
nc.zwave.dimension={'x','y'};

nc=ncshuffle(nc,{'y','x'});

nc=ncstruct(nc);

end

