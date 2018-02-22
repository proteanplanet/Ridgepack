function ridgepack_maptextm(Name,Name2)

% ridgepack_maptextm - This function adds text to a map for a given place name.
%
% function ridgepack_maptextm(Name,Name2)
%
% This function adds text to a map for a given place name. The 
% function looks at an internal database for the name.
%
% Name  - name of place for which a latitude and longitude are
%         required
%
% Name2 - name of second place (optional). If this is provided
%         the function will plot the track between the two 
%         places, Name and Name2, using a great circle.
%
% Andrew Roberts, Naval Postgraduate School, March 2018  (afrobert@nps.edu)
% 

global debug;
if debug; disp(['Entering ',mfilename,'...']); end

if ~ismap(gca)
 error('Must be applied to a current map handle')
end


rot=-90;

[lat,lon]=mplace(Name);

%h=textm(lat,lon,['{\leftarrow}',Name],'Rotation',rot);
h=textm(lat,lon,[' ',Name],'Rotation',rot);
rotatetext(h);
plotm(lat,lon,'r.');

if not(isempty(Name2));

	[lat2,lon2]=mplace(Name2);
	%h=textm(lat2,lon2,['{\leftarrow}',Name2],'Rotation',rot);
	h=textm(lat2,lon2,[' ',Name2],'Rotation',rot);
	rotatetext(h);
	plotm(lat2,lon2,'r.');

	[lattrkgc,lontrkgc] = track2(lat,lon,lat2,lon2);
	plotm(lattrkgc,lontrkgc,'r');

end

drawnow

if debug; disp(['...Leaving ',mfilename]); end

