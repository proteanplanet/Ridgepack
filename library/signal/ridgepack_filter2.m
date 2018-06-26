function [nc,h]=ridgepack_filter2(nc,name,spacedim,Wn,type,n,N,bandwidth)

% ridgepack_filter2 - Classic 2D non-causal FIR filter
%
% function [nc,h]=ridgepack_filter2(nc,name,spacedim,Wn,type,n,N,bandwidth)
% 
% This function filters data in two spatial dimensions for the variable
% 'name' in an nc structure using a hamming-windowed linear-phase FIR 
% non-causal digital filter.  The dimensions over which the filtering is 
% to be done are listed as a cell array in spacedim, and the cutoff
% wavelengths and filter type (low pass, band pass etc) are provided in
% Wn and type, respectively.  The data being filtered may have many dimensions, 
% but only the spacedim dimensions will be filtered, and these must have 
% identical units. The filter order, n, can be specified, or, if n is left empty
% this function finds the minimum order,n, meeting the design criteria set by 
% the filter bandwidth, which is the transition band between passed and stopped 
% frequencies. The stop band attenuation is preset at 30dB (0.1% leakage or less 
% in stop bands) and the pass band ripple is set at 99% or more of the input signal. 
% If the bandwidth is omitted, the transition bands are kept within 10 times
% the frequency resolution of the shortest axis of the two dimensions being 
% filtered.
%
% Note again, as already mentioned, that cutoffs for the 2D filter are specified
% in terms of wavelengths rather then frequencies, which differs from the 
% ridgepack_filter1 function used for time filtering.
%
% If no output arguments are specified, this function acts as a tool for filter
% design, and plots the equivalent frequency and impulse response functions using 
% ncfreqzplot, and also plots the 2D frequency and impulse responses in a second
% figure window. Alternatively, if an output argument is specified, the output 
% is a similar netcdf structure to the input, but with the variable 'name' 
% filtered, and with all variables including dimensions being filtered reduced
% in size by half the filter order (n/2) at each boundary.  As a result, data in 
% the input nc structure should incorporate boundary conditions extending by n/2
% past each boundary if no data loss is expected for the main domain.
%
% Because the filter has been carefully designed to distinguish pass, transition 
% and stop bands, there is no recoloring applied to the pass bands.  However
% this may be required if one chooses to ignore the present passband ripple
% criteria, because this could result in the power of the passband being 
% significantly increased or decreased above or below the input signal.
%
% INPUT:
%
% nc        - nc structure to be filtered
% 
% name      - name of variable in nc structure to be filtered.
%
% spacedim  - cell array of spatial dimensions to be filtered from the data
%             The cell array must have a length of exactly two, being the
%             first and second dimensions of the variable being filtered.
%
% Wn        - Cutoff frequency w1, or a vector giving the stop or pass
%             band filters [w1 w2 w3 ... wn]. The units of the periods
%             must be in the units of the spacedim dimensions.
%
% type      - 'low','high','stop' or 'pass' pass filters.  The default
%             is a low pass filter if this argument is omitted. For
%             the 'stop' and 'pass' band filters, Wn should be a vector
%             [w1 w2 w3 w4 w5] with the first band w1<w<Inf being stopped
%             or passed, respectively, and the second band w2<w<w1 being
%             passed or stopped, respectively.
%
% n         - filter order. If this is left empty [] or omitted, the
%             function will find the optimal order for the given
%             specifications.
%
% N         - maximum number of data points to be used for calculating 
%             the frequency resolution in the grid if less than that
%             supplied along the shortest spatial axis in nc. If this 
%             is left empty [] then N is calculated.
%
% bandwidth - filter bandwidth (ws-wp) where ws is the stop frequency and 
%             wp is the pass frequency for any transition with cutoff 
%             wc=(wp+ws)/2. Default is 10 times the frequency resolution 
%             of the shortest axis of spacedim.  This value should be
%             entered as a multiple of the frequency resolution.
%
%
% OUTPUT:
%
% nc        - nc structure with filtered data for name
%             through each slice of the dimensions spacedime
%
% h         - the filter impulse response function used in filtering.
%
% Ridgepack Version 1.0
% Andrew Roberts, Naval Postgraduate School, March 2018 (afrobert@nps.edu)
% Naval Postgraduate School, Department of Oceanography
%

% check inputs
if nargin<4
  error('Must provide nc structure, name, spacial dimensions and cuttoff wavelength');
elseif not(isstruct(nc));
  error([inputname(1),' is not a structure']);
elseif length(spacedim)~=2
  error('Must specify two dimensions over which filtering is to be done');
end

if nargin<5; type='low'; end

if strcmp(type,'stop') & length(Wn)==1
	error('must specify at least two frequencies for stop band filter');
elseif strcmp(type,'stop') & length(Wn)==2
	type='stop';
elseif strcmp(type,'stop') 
	type='DC-0';
elseif strcmp(type,'pass') & length(Wn)<3
	error('Must specify at least three frequencies for pass band')
elseif strcmp(type,'pass')
	type='DC-1';
elseif ~strcmp(type,{'low','high'})
	error('type must be set to one of: low, high, stop, pass');
end

% get initial dimension order of nc structure
olddimorder=nc.(name).dimension;

% find the number of data points in spatial dimensions by shuffling them last
nc=ncshuffle(nc,spacedim); 
if nargin<7; N=[]; end
if isempty(N)
 N=size(nc.(name).data); 
 N=min(N(end-(length(spacedim)-1):end));
end

% Get frequency axis information for each direction of the grid, but
% use the diagonal resolution which represents the minimum resolution
% over which a wave can be fully resolved in any direction on the grid.
resolution1=abs(nc.(char(spacedim{1})).data(2)-nc.(char(spacedim{1})).data(1));
resolution2=abs(nc.(char(spacedim{2})).data(2)-nc.(char(spacedim{2})).data(1));
fs=1/resolution1; df=fs/N; xunits=nc.(char(spacedim{1})).units;

% check resolution and units are consistent
if ~strcmp(xunits,nc.(char(spacedim{2})).units)
 error('Units of two spatial dimensions are different')
elseif resolution2~=resolution1
 warning(['Resolution in second dimension differs.  It is ',num2str(resolution2),xunits]);
end

disp([char(spacedim{1}),' resolution is ',num2str(resolution1),xunits,...  
      ' with Nyquist frequency ',num2str(fs/2),'/',xunits]);

% set tolerances
if nargin<6; n=[]; end
if nargin<8; 
	bandwidth=10*df; 
else
	bandwidth=bandwidth*df; 
end
attenuation=ridgepack_db(0.010);
ripple=ridgepack_db(0.99);

disp('-------------------------------');
disp('FIR filter');
disp(['Bandwidth tolerance is:   ',num2str(bandwidth)]);
disp(['Attenuation tolerance is: ',num2str(ridgepack_invdb(attenuation))]);
disp(['Ripple tolerance is:      ',num2str(ridgepack_invdb(ripple))]);
disp('-------------------------------');

% convert cutoff wavelengths to frequencies and check against tolerances
Wn=sort(1./Wn);
if any(Wn>fs/2)
	error(['nth elements of Wn greater than the Nyquist frequency: ',...
	       num2str(find(Wn>fs/2))])
elseif bandwidth>=fs/2
	error('bandwidth greater than Nyquist frequency')
end

% determine the order search range, or else set n to the specified even order
if isempty(n)
 nmin=20; nmax=2000;
else
 n=2*ceil(n/2); nmin=n; nmax=n;
end

% search for the best 1D filter along on dimensional path (keeping n even)
for n=nmin:2:nmax;

 % Get difference equation coefficients
 a=1; b=fir1(n,Wn*2/fs,type,hamming(n+1),'scale');

 % get frequency response
 [H,W]=freqz(b,a,n); W=(fs/2)*(W/pi);

 % determine if this frequency response passes criteria
 [pass,description]=ridgepack_fcriteria(H,W,Wn,bandwidth,attenuation,ripple,xunits);

 % if criteria are met, exit from the loop
 if pass; break; end

 % if nmax is reached, exit with an error
 if n==nmax; 
  error(['Filter did not meet specifications up to order ',num2str(n)]); 
 end
  
end 

% convert this 1D information into a 2D impulse and frequency response functions
h =ftrans2(b); 
H2=freqz2(h,n,n); 

% add information to the description
description{1}=['FIR 1D ',char(description{1})];

% either filter the data or provide details of the filter design
if nargout>0

	nc.(name).long_name=['Spatially filtered ',nc.(name).long_name];
	si1=size(nc.(name).data);
	if length(si1)>2; fd=prod(si1(1:end-2)); else; fd=1; end

	nc.(name).data=reshape(nc.(name).data,[fd si1(end-1) si1(end)]);

	for i=1:fd;
         disp(['Filtering record ',num2str(i)])
	 nc.(name).data(i,:,:)=filter2(h,squeeze(nc.(name).data(i,:,:)),'same');
        end

	nc.(name).data=reshape(nc.(name).data,si1);

	[variablenames,numbervariables]=ncname(nc);
	for i=1:2
         dimname=char(spacedim{i});
         nc=ncshuffle(nc,{dimname});
	 for m = 1:numbervariables
	   name=char(variablenames(m));
	   if ~isempty(intersect({dimname},nc.(name).dimension)) & isfield(nc.(name),'data')
	    si1=size(nc.(name).data);
	    if si1(2)==1; si1=si1(end:-1:1); end
	    nc.(name).data=reshape(nc.(name).data,[prod(si1(1:end-1)),si1(end)]);
	    nc.(name).data=nc.(name).data(:,1+n/2:end-n/2);
	    si2=size(nc.(name).data);
	    nc.(name).data=reshape(nc.(name).data,[si1(1:end-1) si2(end)]);
           end
	 end
        end

	if isfield(nc.attributes,'comment')
	  nc.attributes.comment=[nc.attributes.comment,' Spatially filtered using ',...
	                         nccellcat(description)];
        else
	  nc.attributes.comment=['Spatially filtered using ',nccellcat(description)];
        end

	nc=ncshuffle(nc,olddimorder);

	nc=ncstruct(nc);

else
	% plot filter characteristics on a graph
	figure(1); clf; 
	ncfreqzplot(H,W,b,a,description,xunits);

	% plot 2D filter characteristics on a mesh
	figure(2); clf; 
	[H2,W1,W2]=freqz2(h,n,n); W1=(fs/2)*W1; W2=(fs/2)*W2;

	subplot(211), mesh(W1,W2,H2);
	xlim([min(W1) max(W1)]); ylim([min(W2) max(W2)]); zlim([min(H2(:)) max(H2(:))]); 
	xlabel(['w_x (',xunits,'^{-1})']); ylabel(['w_y (',xunits,'^{-1})'])
	title('Frequency response H(w_x,w_y)')

	subplot(212), mesh([1:1:size(h,1)],[1:1:size(h,2)],h); 
	xlim([1 size(h,1)]); ylim([1 size(h,2)]); zlim([min(h(:)) max(h(:))]); 
	xlabel('x'); ylabel('y')
	title('Impulse response h(x,y)')

	clear nc;
end

disp('-------------------------------');
disp(['Maximum frequency response: ',num2str(max(abs(H2(:))))])
disp(['Minimum frequency response: ',num2str(min(abs(H2(:))))])
disp('-------------------------------');

end % function
