function [pass,description]=ridgepack_fcriteria(H,W,Wn,bandwidth,attenuation,ripple,units)

% ridgepack_fcriteria - Determines if a given filter exceeds given design criteria
%
% function [pass,description]=ridgepack_fcriteria(H,W,Wn,bandwidth,attenuation,ripple)
%
% This function determines whether the nth-order digital filter with 
% frequency response H and discrete frequencies W in cycles per input units
% fits the desired criteria of cutoff frequencies Wn, transition bandwidth, 
% stop band attenuation and passband ripple.
% 
% INPUT:
%
% H           - Frequency response
%
% W           - Discrete frequency axis in cycles/units
%
% Wn          - Vector of cuttoff frequencies defining 
%               the transition between pass and stop bands
%               This is the same input as required for FIR1.
%
% bandwidth   - bandwidth in cycles per units
%
% attenuation - stop band attenuation in db
%
% ripple      - passband ripple in db
%
% units       - time or space units being filtered (days, months, km etc)
%
%
% OUTPUT:
%
% pass        - logical stating whether or not the provided frequency
%               response passes the given criteria.
%
% description - string stating why H failed the test, or else
%               listing the characteristics of a filter that
%               passes the test.
%
% Ridgepack Version 1.0
% Andrew Roberts, Naval Postgraduate School, March 2018 (afrobert@nps.edu)
% Naval Postgraduate School, Department of Oceanography
%

% Find the index values of the pass and stop bands 
passi=find(abs(H)>=ridgepack_invdb(ripple));
stopi=find(abs(H)<=ridgepack_invdb(attenuation));

% Make sure upper ripple bound is not exceeded
if ~isempty(find(abs(H)>2-ridgepack_invdb(ripple)));
	description='Upper ripple bound exceeded';
	pass=false; return;
end

% check that criteria have been reached
if isempty(passi) | isempty(stopi)
	description='Pass or stop bands do not exist';
	pass=false; return;
end

% Get the actual minimum stopband attenuation and actual passband ripple 
attenuation=-ridgepack_db(abs(max(abs(H(stopi)))-min(abs(H(stopi))))); 
ripple=-ridgepack_db(abs(min(abs(H(passi))))); 

% Determine the frequency of each stop band using the index values
Wstop(1)=W(stopi(1));
k=1;
for i=2:length(stopi)
	if stopi(i)~=stopi(i-1)+1
		k=k+1; Wstop(k)=W(stopi(i-1));
		k=k+1; Wstop(k)=W(stopi(i));
	end
end
Wstop(length(Wstop)+1)=W(stopi(end));

% Determine the frequency of each pass band using the index values
Wpass(1)=W(passi(1));
k=1;
for i=2:length(passi)
	if passi(i)~=passi(i-1)+1
		k=k+1; Wpass(k)=W(passi(i-1));
		k=k+1; Wpass(k)=W(passi(i));
	end
end
Wpass(length(Wpass)+1)=W(passi(end));

% Determine the number of transistion bands in the real filter
passn=length(find(Wpass>W(1) & Wpass<W(end)));
stopn=length(find(Wstop>W(1) & Wstop<W(end)));

% Check to see if pass and stop bands satisfy fundamental criteria
if W(passi(end)) ~= W(end) & W(stopi(end)) ~= W(end)
	description='Pass and stop bands do not cover full frequency range';
	pass=false; return; 
elseif passn~=stopn
	description='Incorrect number of pass or stop bands';
	pass=false; return; 
elseif passn~=length(Wn)
	description='Number of transition bands different from Wn';
	pass=false; return; 
end

% Generate transition bandwidths and text information on the 
% pass and stop bands in both frequency and time domains
if Wpass(1)<Wstop(1)
	for i=1:2:100000;
		if i==100000 ; error('transband overflow') ; end
		if length(Wpass)<i+1 | length(Wstop)<i ; break; end
		if i>length(Wn)
			description='Number of bands exceeds those in Wn';
			pass=false;
			return;
		elseif Wn(i)>Wstop(i) | Wn(i)<Wpass(i+1)
			description='Wn outside of specified band';
			pass=false;
			return;
		end
		transband(i)=Wstop(i)-Wpass(i+1);
		if length(Wpass)<i+2 | length(Wstop)<i+1 ; break; end
		if i+1>length(Wn)
			description='Number of bands exceeds those in Wn';
			pass=false;
			return;
		elseif Wn(i+1)<Wstop(i+1) | Wn(i+1)>Wpass(i+2)
			description='Wn outside of specified band';
			pass=false;
			return; 
		end
		transband(i+1)=Wpass(i+2)-Wstop(i+1);
	end
elseif Wpass(1)>Wstop(1)
	for i=1:2:100000;
		if i==100000 ; error('transband overflow') ; end
		if length(Wstop)<i+1 | length(Wpass)<i ; break; end
		if i>length(Wn)
			description='Number of bands exceeds those in Wn';
			pass=false;
			return; 
		elseif Wn(i)<Wstop(i+1) | Wn(i)>Wpass(i)
			description='Wn outside of specified band';
			pass=false;
			return;
		end
		transband(i)=Wpass(i)-Wstop(i+1);
		if length(Wstop)<i+2 | length(Wpass)<i+1 ; break; end
		if i+1>length(Wn)
			description='Number of bands exceeds those in Wn';
			pass=false;
			return; 
		elseif Wn(i+1)>Wstop(i+2) | Wn(i+1)<Wpass(i+1)
			description='Wn outside of specified band';
			pass=false;
			return;
		end
		transband(i+1)=Wstop(i+2)-Wpass(i+1);
	end
else
	description='Pass and stop bands are coincident';
	pass=false; return; 
end

% check that that the bandwidth is within the desired specification
if ~all(transband<=bandwidth)
	description='Bandwidth exceeds criteria';
	pass=false; return; 
end

% if this point has been reached, the filter has passed the test
pass=true;

% generate bandwidth information since the filter has now passed all criteria
if Wpass(1)<Wstop(1)
	bandsf=['Bands (cycles/',units,'): '];
	for i=1:2:100000;
		if i==100000 ; error('band overflow') ; end
		if length(Wpass)<i+1 | length(Wstop)<i ; break; end
		if i>length(Wpass); break; end
		bandsf=[bandsf,num2str(Wpass(i)),'|pass|',num2str(Wpass(i+1)),'  '];
		if i>length(Wstop); break; end
		bandsf=[bandsf,num2str(Wstop(i)),'|stop|',num2str(Wstop(i+1)),'  '];
	end
	if Wpass(2)>1.0 & strcmp(units,'days');
		Wpass=24./Wpass;
		Wstop=24./Wstop;
		bandsp=['Bands (hours): '];
	else
		Wpass=1./Wpass;
		Wstop=1./Wstop;
		bandsp=['Bands (',units,'): '];
	end
	for i=1:2:100000;
		if i==100000 ; error('band overflow') ; end
		if i>length(Wpass); break; end
		bandsp=[bandsp,num2str(Wpass(i)),'|pass|',num2str(Wpass(i+1)),'  '];
		if i>length(Wstop); break; end
		bandsp=[bandsp,num2str(Wstop(i)),'|stop|',num2str(Wstop(i+1)),'  '];
	end

elseif Wpass(1)>Wstop(1)
	bandsf=['Bands (cycles/',units,'): '];
	for i=1:2:100000;
		if i==100000 ; error('band overflow') ; end
		if i>length(Wstop); break; end
		bandsf=[bandsf,num2str(Wstop(i)),'|stop|',num2str(Wstop(i+1)),'  '];
		if i>length(Wpass); break; end
		bandsf=[bandsf,num2str(Wpass(i)),'|pass|',num2str(Wpass(i+1)),'  '];
	end
	if Wpass(2)>1.0 & strcmp(units,'days');
		Wpass=24./Wpass;
		Wstop=24./Wstop;
		bandsp=['Bands (hours): '];
	else
		Wpass=1./Wpass;
		Wstop=1./Wstop;
		bandsp=['Bands (',units,'): '];
	end
	for i=1:2:100000;
		if i==100000 ; error('band overflow') ; end
		if i>length(Wstop); break; end
		bandsp=[bandsp,num2str(Wstop(i)),'|stop|',num2str(Wstop(i+1)),'  '];
		if i>length(Wpass); break; end
		bandsp=[bandsp,num2str(Wpass(i)),'|pass|',num2str(Wpass(i+1)),'  '];
	end
end

% Add this information to the description
info=['Filter order=',num2str(length(W)),...
      ', Passband Ripple=',num2str(ripple),...
      'db, Stopband Attenuation=',num2str(attenuation),'db'];


% Provide information on the filter 
description={info,bandsf,bandsp};
disp(info);
disp(bandsf);
disp(bandsp);


end % function
