function [xlats,xlons,xverts]=...
                   ridgepack_e3smcontour(tlats,tlons,tverts,cverts)


% Find NaN's seperating the vertices. idx will have a length
% of n-1 where n is the number of closed loops.
tvertsnan=[NaN tverts NaN];
idx=find(isnan(tvertsnan));

% Also find the number of coastal closed loops
cvertsnan=[NaN cverts NaN];
cdx=find(isnan(cvertsnan));

% create new set of contours 
xvertsnan=tvertsnan;
xlatsnan=[NaN tlats NaN];
xlonsnan=[NaN tlons NaN];

nbreaks=0; % number of line breaks that need to be inserted
breakvert=[]; % break of vertex index

for i=1:length(idx)-1

 % vertices along a close loop
 vertcircleidx=[idx(i)+1:idx(i+1)-1 idx(i)+2:idx(i+1)-1];
 vertcircle=tvertsnan(vertcircleidx);

 for k=1:length(cdx)-1

  coastcirclecdx=[cdx(k)+1:cdx(k+1)-1 cdx(k)+2:cdx(k)+3];
  coastcircle=cvertsnan(coastcirclecdx);

  for j=1:length(vertcircle)-2

   if ~isempty(strfind(coastcircle,vertcircle(j:j+1))) | ...
      ~isempty(strfind(coastcircle,vertcircle(j+1:-1:j)))
     nbreaks=nbreaks+1;
     breakvert(nbreaks)=vertcircleidx(j);
   end

  end

 end
end

% Now scan for shared sides between the new contours and 
% the coastline, and place a NaN in between
sortidx=sort(breakvert);
for i=length(sortidx):-1:1
 xvertsnan=[xvertsnan(1:sortidx(i)) NaN xvertsnan(sortidx(i)+1:end)];
 xlatsnan=[xlatsnan(1:sortidx(i)) NaN xlatsnan(sortidx(i)+1:end)];
 xlonsnan=[xlonsnan(1:sortidx(i)) NaN xlonsnan(sortidx(i)+1:end)];
end

% Now remove sequences of [NaN Number NaN] and replace with NaN
% Also removing the NaN that padded the start and finish of the 
% sequence in the previous search, and sequences of NaNs.

xverts=xvertsnan;
xlats=xlatsnan;
xlons=xlonsnan;

% First, fill in single vertices with NaNs
for i=2:length(xverts)-1
 if (isnan(xverts(i-1)) & ...
    ~isnan(xverts(i)) & ...
     isnan(xverts(i+1)))
  xvertsnan(i)=NaN;
  xlatsnan(i)=NaN;
  xlonsnan(i)=NaN;
 end
end

xverts=[];
xlats=[];
xlons=[];

% Now, replace repeated NaNs with a single one
k=0;
for i=2:length(xvertsnan)
 if ~(isnan(xvertsnan(i-1)) & isnan(xvertsnan(i)))
  k=k+1;
  xverts(k)=xvertsnan(i);
  xlats(k)=xlatsnan(i);
  xlons(k)=xlonsnan(i);
 end
end

% quick error check
for i=2:length(xverts)
 if (isnan(xverts(i-1)) & isnan(xverts(i)))
  error('NaNs repeated')
 end
end




