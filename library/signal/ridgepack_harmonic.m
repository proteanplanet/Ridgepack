function [nc]=ridgepack_harmonic(nc,name,lengtht,constituent)

% ridgepack_harmonic - Calculate the phase and amplitude of a field with constant frequency
%
% function [nc]=ridgepack_harmonic(nc,name,lengtht,constituent)
%
% This function calculates the phase and amplitude of a 2D sinusoidal waveform 
% (2 spatial dimensions and 1 time dimension) with constant frequency in the input 
% nc structure. Phase may be referenced to a particular tidal constituent if required.
% 
% INPUT:
%
% nc - nc structure with 'name' as a variable
% name - variable in the nc structure with a 2D constant frequency waveform
% lengtht - length, in number of samples, of the waveform.
% constituent - optional string giving a tidal component against which the 
%               phase may be referenced. eg. constituent=['M2']; 
%	        This provides an approximation only if used with tidal magnitude.
%
%
% OUTPUT:
%
% nc - nc structure with the phase and amplitude of name
%
% Ridgepack Version 1.0
% Andrew Roberts, Naval Postgraduate School, March 2018 (afrobert@nps.edu)
% Naval Postgraduate School, Department of Oceanography
%

% input check
if ~isstruct(nc)
 error('nc is not a structure')
elseif ~isfield(nc,'time')
 error('time is not a variable of nc')
elseif ~isfield(nc,name)
 error([name,' is not a variable of nc'])
elseif length(nc.time.data)<round(lengtht/4)
 error('You must provide at least a quarter of a period of the waveform')
end

% fill in nc structure details
nc.([name,'_phase'])=nc.(name);
nc.([name,'_phase']).long_name=[name,' phase'];
nc.([name,'_phase']).units='degrees';
nc.([name,'_phase']).dimension=ncsetdiff(nc.(name).dimension,{'time'});

% fill in nc structure details
nc.([name,'_amplitude'])=nc.(name);
nc.([name,'_amplitude']).long_name=[name,' amplitude'];
nc.([name,'_amplitude']).dimension=nc.([name,'_phase']).dimension;

% place time as last dimension and get number of last (time) dimension
nc=ncshuffle(nc,{'time'});
nd=ndims(nc.(name).data);

% get the nearest time index value to the reference time
if nargin<4 | isempty(constituent)

 phasechange=0;

elseif strcmpi(constituent,'M2')

 % m2 period in days
 m2period=12.4206/24;
 disp(['     M2 period in days is: ',num2str(m2period)])

 % zero phase time for greenwhich
 fiducial=datenum([1996,7,20,0,26,40]);

 % zero phase within time period close to the start of supplied record
 fiducial=fiducial+m2period*ceil((nc.time.data(1)-fiducial)/m2period);

 % display information
 disp(['Nearest fiducial point is: ',datestr(fiducial)])
 disp(['            Start time is: ',datestr(nc.time.data(1))])
 disp(['              End time is: ',datestr(nc.time.data(end))])

 % error checking
 if fiducial<nc.time.data(1) ; error('Fiducial calculation error 1'); end

 % get phase subtraction to be made to each calculation 
 phasechange=360*((nc.time.data(1)-fiducial)/m2period);

 % zero phase correction
 disp(['  Zero phase corrected to: ',num2str(phasechange),' degrees'])

else

 error('Specified constituent not recognized')

end

% get the index at start time and a quarter period later
t0=1; tq=t0+round(lengtht/4);
disp(['Quarter cycle occurs at record ',num2str(tq)]);

% calculate phase and amplitude
nc.([name,'_phase']).data=mod(wrapTo360(rad2deg(atan2(nc.(name).data(:,:,tq),nc.(name).data(:,:,t0))))...
                          + phasechange,360);
nc.([name,'_amplitude']).data=sqrt(nc.(name).data(:,:,tq).^2 + nc.(name).data(:,:,t0).^2);

[nc,out]=ncstruct(nc);


