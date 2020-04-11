function [nc]=ridgepack_pthresholde3sm(ncvert,ncc,varc,threshold,...
                                       centlat,centlon,horizon)

% ridgepack_pthresholde3sm - generate a threshold on an unstructured mesh
%
% function [nc]=ridgepack_pthresholde3sm(ncc,varc,threshold,ncvert,centlat,centlon,horizon)
%
% This function .... description
%  
% INPUT:
%
% ncc       - netcdf structure including
% varc      - variable in ncc to be used
% threshold - threshold value used as a boundary on unstructured mesh
% ncvert    - vertices on the unstructed mesh
% centlat   - center latitude of plotted satellite view (optional)
% centlon   - center longitude if plotted in satellite view (optional)
% horizon   - horizon, in degrees of satellite view (optional)
%
% Note that centlat, centlon and horizon should not be provided
% if nc is specificied.
%
% 
% OUTPUT:
% 
% nc - netcdf structure with lats, longs, mesh cell and vertex
%      indices
%
% Ridgepack Version 1.1
% Andrew Roberts, Los Alamos National Laboratory, April 3, 2020 (afroberts@lanl.gov)



if nargin<4
 error('There is no threshold value set, and possibly more')
end

if nargin<5
 centlat=90;
elseif nargin>3 & nargout==1
 error('you are asking to both plot and produce coastal data')
elseif ~isnumeric(centlat)
 error('centlat should be a real number between -90 and 90')
elseif centlat>90 | centlat<-90
 error('centlat should be between -90 and 90')
end

if nargin<6
 centlon=0;
elseif ~isnumeric(centlon)
 error('centlon should be a real number between -90 and 90')
elseif centlat>180 | centlat<-180
 error('centlon should be between -180 and 180')
end

if nargin<7
 horizon=90;
elseif ~isnumeric(horizon)
 error('horizon must be a real number between 0 and 90')
elseif horizon<0 | horizon>90
 error('horizon must be between 0 and 90')
end

% get indices of cells that match the criteria
cidx=find(ncc.(varc).data>threshold);
idx=[1:length(ncvert.nCells.data)]';

% get size of total vertices in cells within criteria
cnlength=sum(ncvert.nEdgesOnCell.data(cidx));
nlength=sum(ncvert.nEdgesOnCell.data(:));

% now use it to set the of vertexlist
cvertexlist=zeros([cnlength 1]);
vertexlist=zeros([nlength 1]);

% compile list of vertices
k=1;
for i=cidx'
 cvertexlist(k:k+ncvert.nEdgesOnCell.data(i)-1)=...
  ncvert.verticesOnCell.data(1:ncvert.nEdgesOnCell.data(i),i);
 k=k+ncvert.nEdgesOnCell.data(i);
end

k=1;
for i=idx'
 vertexlist(k:k+ncvert.nEdgesOnCell.data(i)-1)=...
  ncvert.verticesOnCell.data(1:ncvert.nEdgesOnCell.data(i),i);
 k=k+ncvert.nEdgesOnCell.data(i);
end

% sort list of data to find number of occurrences
csortexvertex=sort(cvertexlist);
sortexvertex=sort(vertexlist);

% find number of c cells per vertex
disp('Doing length search')
cdixv=zeros([length(ncvert.nVertices.data) 1]);
for i=1:length(csortexvertex)
 cdixv(csortexvertex(i))=cdixv(csortexvertex(i))+1;
end
dixv=zeros([length(ncvert.nVertices.data) 1]);
for i=1:length(sortexvertex)
 dixv(sortexvertex(i))=dixv(sortexvertex(i))+1;
end

% now calculate mask of vertices at land border noting that 
% a only a vertex bording 1 or 2 cells sits on the coast.
cpassx=NaN*zeros([length(ncvert.nVertices.data) 1]);
cpassx(cdixv==1 | cdixv==2)=1;
passx=NaN*zeros([length(ncvert.nVertices.data) 1]);
passx(dixv==1 | dixv==2)=1;

% now find where cells have vertex bordering 2 cells
cpassy=NaN*zeros([length(ncvert.nVertices.data) 1]);
cpassy(cdixv==2)=1;
passy=NaN*zeros([length(ncvert.nVertices.data) 1]);
passy(dixv==2)=1;

% mask out lats and longs that are not at threshold edge
ncvert.clatitude=ncvert.latitude;
ncvert.clongitude=ncvert.longitude;
ncvert.clatitude.data=cpassx.*ncvert.clatitude.data;
ncvert.clongitude.data=cpassx.*ncvert.clongitude.data;

% mask out areas not in the threshold plotting area
if nargout==0
 disp('Reducing to just the threshold plotting area')
 maxth=deg2rad(horizon);
 for i=cidx'

  maxidx=ncvert.nEdgesOnCell.data(i);

  la=ncvert.clatitude.data(ncvert.verticesOnCell.data(1:maxidx,i));
  lo=ncvert.clongitude.data(ncvert.verticesOnCell.data(1:maxidx,i));

  [x,y,z,ph,th]=ridgepack_satfwd(rad2deg(la),rad2deg(lo),...
                                centlat,centlon,horizon,1);

  % filter cells not in frame, and find cropping limit
  if all(isnan(x))
   ncvert.clatitude.data(ncvert.verticesOnCell.data(1:maxidx,i))=NaN;
   ncvert.clongitude.data(ncvert.verticesOnCell.data(1:maxidx,i))=NaN;
  end

 end
end

% calculate mask for cells on the threshold edge and for 
% sides with 2 adjoined vertices. Flagcell indicates
% a cell in which the cell has two adjacent 2-vertex
% locations, and flags those vertices for a later search. 
% This is to stop them being joined on the coastline.
disp('Masking to threshold edge cells only')
cpasscell=zeros([length(ncvert.nCells.data(cidx)) 1]);
k=1;
for i=cidx'
 if any(cpassx(ncvert.verticesOnCell.data...
       (1:ncvert.nEdgesOnCell.data(i),i))==1) 
  cpasscell(i)=1;
 end
 for j=1:ncvert.nEdgesOnCell.data(i)
  if j<ncvert.nEdgesOnCell.data(i) & ...
     all(cpassy(ncvert.verticesOnCell.data(j:j+1,i))==1) 
   cflagcell(k)=i;
   vflag=sort(ncvert.verticesOnCell.data(j:j+1,i));
   cvertexflag1(k)=vflag(1);
   cvertexflag2(k)=vflag(2);
   k=k+1;
  elseif j==ncvert.nEdgesOnCell.data(i) & ...
     all(cpassy(ncvert.verticesOnCell.data([1 j],i))==1)
   cflagcell(k)=i;
   vflag=sort(ncvert.verticesOnCell.data([1 j],i));
   cvertexflag1(k)=vflag(1);
   cvertexflag2(k)=vflag(2);
   k=k+1;
  end
 end
end
cidxn=find(cpasscell==1);

disp('Masking to threshold edge cells only')
passcell=zeros([length(ncvert.nCells.data(idx)) 1]);
k=1;
for i=idx'
 if any(passx(ncvert.verticesOnCell.data...
       (1:ncvert.nEdgesOnCell.data(i),i))==1)
  passcell(i)=1;
 end
 for j=1:ncvert.nEdgesOnCell.data(i)
  if j<ncvert.nEdgesOnCell.data(i) & ...
     all(passy(ncvert.verticesOnCell.data(j:j+1,i))==1)
   flagcell(k)=i;
   flag=sort(ncvert.verticesOnCell.data(j:j+1,i));
   vertexflag1(k)=flag(1);
   vertexflag2(k)=flag(2);
   k=k+1;
  elseif j==ncvert.nEdgesOnCell.data(i) & ...
     all(passy(ncvert.verticesOnCell.data([1 j],i))==1)
   flagcell(k)=i;
   flag=sort(ncvert.verticesOnCell.data([1 j],i));
   vertexflag1(k)=flag(1);
   vertexflag2(k)=flag(2);
   k=k+1;
  end
 end
end
idxn=find(passcell==1);

% find shared edges on cells of ice edge
[cvflag1,I]=sort(cvertexflag1);
cvflag2=cvertexflag2(I);
cfcell=cflagcell(I);

k=1;
for i=1:length(cvflag1)
 csearchj=false;
 for j=1:length(cvflag1)
  if cvflag1(i)==cvflag1(j) & cvflag2(i)==cvflag2(j) & i~=j 
   cflagcoast1(k)=cvflag1(i);
   cflagcoast2(k)=cvflag2(i);
   ccellflag1(k)=cfcell(i);
   ccellflag2(k)=cfcell(j);
   k=k+1;
   csearchj=true;
  end
  if csearchj
   break
  end
 end
end

% find shared edges on cells of ice edge
[vflag1,I]=sort(vertexflag1);
vflag2=vertexflag2(I);
fcell=flagcell(I);

k=1;
for i=1:length(vflag1)
 searchj=false;
 for j=1:length(vflag1)
  if vflag1(i)==vflag1(j) & vflag2(i)==vflag2(j) & i~=j 
   flagcoast1(k)=vflag1(i);
   flagcoast2(k)=vflag2(i);
   cellflag1(k)=fcell(i);
   cellflag2(k)=fcell(j);
   k=k+1;
   searchj=true;
  end
  if searchj
   break
  end
 end
end

% weed out coastal sides of cells
for k=1:length(cflagcoast1);
 for j=1:length(flagcoast1);
  if (cflagcoast1(k)==flagcoast1(j) & cflagcoast2(k)==flagcoast2(j))
   xflagcoast(k)=NaN;
   xcellflag(k)=NaN;
  else 
   xflagcoast(k)=1;
   xcellflag(k)=1;
  end
 end
end


% only pass through threshold cells
disp('Calculating the extent line')
if length(cidxn)>0

      px0=[];
      py0=[];
      pcells0=[];
      pverts0=[];

      px1=[];
      py1=[];
      pcells1=[];
      pverts1=[];

      px2=[];
      py2=[];
      pcells2=[];
      pverts2=[];

      px3=[];
      py3=[];
      pcells3=[];
      pverts3=[];

      for i=1:length(cidxn)

       maxidx=ncvert.nEdgesOnCell.data(cidxn(i));

       if (any(ccellflag1==cidxn(i)) | any(ccellflag2==cidxn(i))) 

        if any(ccellflag1==cidxn(i))
         k=find(ccellflag1==cidxn(i));
        else
         k=find(ccellflag2==cidxn(i));
        end

        if size(k,2)==4 

         error('Encountered four shared coastal sides')

        elseif size(k,2)==3 

         idx1=find(ncvert.verticesOnCell.data(1:maxidx,cidxn(i))==...
                  cflagcoast1(k(1)))*xflagcoast(k(1));

         idx2=find(ncvert.verticesOnCell.data(1:maxidx,cidxn(i))==...
                  cflagcoast2(k(1)))*xflagcoast(k(1));

         idx3=find(ncvert.verticesOnCell.data(1:maxidx,cidxn(i))==...
                  cflagcoast1(k(2)))*xflagcoast(k(2));

         idx4=find(ncvert.verticesOnCell.data(1:maxidx,cidxn(i))==...
                  cflagcoast2(k(2)))*xflagcoast(k(2));

         idx5=find(ncvert.verticesOnCell.data(1:maxidx,cidxn(i))==...
                  cflagcoast1(k(3)))*xflagcoast(k(3));

         idx6=find(ncvert.verticesOnCell.data(1:maxidx,cidxn(i))==...
                  cflagcoast2(k(3)))*xflagcoast(k(3));

         xxp1=sort([idx1 idx2]);
         xxp2=sort([idx3 idx4]);
         xxp3=sort([idx5 idx6]);

         if xxp1(1)==1 & xxp1(2)==maxidx 
          idx5=xxp1(1);
          idx6=xxp1(2);
          if xxp2(1)<xxp3(1)
           idx1=xxp2(1);
           idx2=xxp2(2);
           idx3=xxp3(1);
           idx4=xxp3(2);
          else
           idx1=xxp3(1);
           idx2=xxp3(2);
           idx3=xxp2(1);
           idx4=xxp2(2);
          end
         elseif xxp2(1)==1 & xxp2(2)==maxidx
          idx5=xxp2(1);
          idx6=xxp2(2);
          if xxp1(1)<xxp3(1)
           idx1=xxp1(1);
           idx2=xxp1(2);
           idx3=xxp3(1);
           idx4=xxp3(2);
          else
           idx1=xxp3(1);
           idx2=xxp3(2);
           idx3=xxp1(1);
           idx4=xxp1(2);
          end
         elseif xxp3(1)==1 & xxp3(2)==maxidx
          idx5=xxp3(1);
          idx6=xxp3(2);
          if xxp1(1)<xxp2(1)
           idx1=xxp1(1);
           idx2=xxp1(2);
           idx3=xxp2(1);
           idx4=xxp2(2);
          else
           idx1=xxp2(1);
           idx2=xxp2(2);
           idx3=xxp1(1);
           idx4=xxp1(2);
          end
         elseif xxp1(1)<xxp2(1) & xxp2(1)<xxp3(1)
          idx1=xxp1(1);
          idx2=xxp1(2);
          idx3=xxp2(1);
          idx4=xxp2(2);
          idx5=xxp3(1);
          idx6=xxp3(2);
         elseif xxp1(1)<xxp3(1) & xxp3(1)<xxp2(1)
          idx1=xxp1(1);
          idx2=xxp1(2);
          idx3=xxp3(1);
          idx4=xxp3(2);
          idx5=xxp2(1);
          idx6=xxp2(2);
         elseif xxp3(1)<xxp1(1) & xxp1(1)<xxp2(1)
          idx1=xxp3(1);
          idx2=xxp3(2);
          idx3=xxp1(1);
          idx4=xxp1(2);
          idx5=xxp2(1);
          idx6=xxp2(2);
         elseif xxp3(1)<xxp2(1) & xxp2(1)<xxp1(1)
          idx1=xxp3(1);
          idx2=xxp3(2);
          idx3=xxp2(1);
          idx4=xxp2(2);
          idx5=xxp1(1);
          idx6=xxp1(2);
         elseif xxp2(1)<xxp1(1) & xxp1(1)<xxp3(1)
          idx1=xxp2(1);
          idx2=xxp2(2);
          idx3=xxp1(1);
          idx4=xxp1(2);
          idx5=xxp3(1);
          idx6=xxp3(2);
         elseif xxp2(1)<xxp3(1) & xxp3(1)<xxp1(1)
          idx1=xxp2(1);
          idx2=xxp2(2);
          idx3=xxp3(1);
          idx4=xxp3(2);
          idx5=xxp1(1);
          idx6=xxp1(2);
         end

         lat(1:idx1)=ncvert.clatitude.data(...
                ncvert.verticesOnCell.data(...
                1:idx1,cidxn(i)))*180/pi;
         lon(1:idx1)=ncvert.clongitude.data(...
                ncvert.verticesOnCell.data(...
                1:idx1,cidxn(i)))*180/pi;

         verts(1:idx1)=ncvert.verticesOnCell.data(...
                                     1:idx1,cidxn(i));
         cells(1:idx1)=cidxn(1);

         lat(idx1+1)=NaN;
         lon(idx1+1)=NaN;

         verts(idx1+1)=NaN;
         cells(idx1+1)=cidxn(i);

         lat(idx2+1:idx3+1)=ncvert.clatitude.data(...
                ncvert.verticesOnCell.data(...
                idx2:idx3,cidxn(i)))*180/pi;
         lon(idx2+1:idx3+1)=ncvert.clongitude.data(...
                ncvert.verticesOnCell.data(...
                idx2:idx3,cidxn(i)))*180/pi;

         verts(idx2+1:idx3+1)=ncvert.verticesOnCell.data(...
                                 idx2:idx3,cidxn(i));
         cells(idx2+1:idx3+1)=cidxn(i);

         lat(idx3+2)=NaN;
         lon(idx3+2)=NaN;

         verts(idx3+2)=NaN;
         cells(idx3+2)=cidxn(i);

         if idx5==1

          lat(idx4+2:maxidx+2)=ncvert.clatitude.data(...
                ncvert.verticesOnCell.data(...
                idx4:maxidx,cidxn(i)))*180/pi;
          lon(idx4+2:maxidx+2)=ncvert.clongitude.data(...
                ncvert.verticesOnCell.data(...
                idx4:maxidx,cidxn(i)))*180/pi;

          verts(idx4+2:maxidx+2)=ncvert.verticesOnCell.data(...
                                          idx4:maxidx,cidxn(i));
          cells(idx4+2:maxidx+2)=cidxn(i);

         elseif idx5>idx4

          lat(idx4+2:idx5+2)=ncvert.clatitude.data(...
                ncvert.verticesOnCell.data(...
                idx4:idx5,cidxn(i)))*180/pi;
          lon(idx4+2:idx5+2)=ncvert.clongitude.data(...
                ncvert.verticesOnCell.data(...
                idx4:idx5,cidxn(i)))*180/pi;

          verts(idx4+2:idx5+2)=ncvert.verticesOnCell.data(...
                                          idx4:idx5,cidxn(i));
          cells(idx4+2:idx5+2)=cidxn(i);

          lat(idx5+3)=NaN;
          lon(idx5+3)=NaN;

          verts(idx5+3)=NaN;
          cells(idx5+3)=NaN;

          lat(idx6+3:maxidx+3)=ncvert.clatitude.data(...
                ncvert.verticesOnCell.data(...
                idx6:maxidx,cidxn(i)))*180/pi;
          lon(idx6+3:maxidx+3)=ncvert.clongitude.data(...
                ncvert.verticesOnCell.data(...
                idx6:maxidx,cidxn(i)))*180/pi;

          verts(idx6+3:maxidx+3)=ncvert.verticesOnCell.data(...
                                           idx6:maxidx,cidxn(i));
          cells(idx6+3:maxidx+3)=cidxn(i);

          lon(maxidx+4)=lon(1);
          lat(maxidx+4)=lat(1);

          verts(maxidx+4)=verts(1);
          cells(maxidx+4)=cidxn(i);

         end

         % append to string of points
         px0=[px0; NaN; lat(:)];
         py0=[py0; NaN; lon(:)];
         pverts0=[pverts0; NaN; verts(:)];
         pcells0=[pcells0; NaN; cells(:)];

         clear lat lon verts cells idx1 idx2 idx3 idx4 idx5 idx6

        elseif size(k,2)==2 

         idx1=find(ncvert.verticesOnCell.data(1:maxidx,cidxn(i))==...
                  cflagcoast1(k(1)))*xflagcoast(k(1));

         idx2=find(ncvert.verticesOnCell.data(1:maxidx,cidxn(i))==...
                  cflagcoast2(k(1)))*xflagcoast(k(1));

         idx3=find(ncvert.verticesOnCell.data(1:maxidx,cidxn(i))==...
                  cflagcoast1(k(2)))*xflagcoast(k(2));

         idx4=find(ncvert.verticesOnCell.data(1:maxidx,cidxn(i))==...
                  cflagcoast2(k(2)))*xflagcoast(k(2));

         xxp1=sort([idx1 idx2]);
         xxp2=sort([idx3 idx4]);

         idx1=xxp1(1);
         idx2=xxp1(2);
         idx3=xxp2(1);
         idx4=xxp2(2);

         if idx3==1 & idx4==maxidx
          idx1=xxp1(1);
          idx2=xxp1(2);
          idx3=xxp2(1);
          idx4=xxp2(2);
         elseif idx1==1 & idx2==maxidx
          idx1=xxp2(1);
          idx2=xxp2(2);
          idx3=xxp1(1);
          idx4=xxp1(2);
         elseif idx2<idx3
          idx1=xxp1(1);
          idx2=xxp1(2);
          idx3=xxp2(1);
          idx4=xxp2(2);
         elseif idx4<idx1
          idx1=xxp2(1);
          idx2=xxp2(2);
          idx3=xxp1(1);
          idx4=xxp1(2);
         else
          error('Configuration not found')
         end

         lat(1:idx1)=ncvert.clatitude.data(...
                ncvert.verticesOnCell.data(...
                1:idx1,cidxn(i)))*180/pi;
         lon(1:idx1)=ncvert.clongitude.data(...
                ncvert.verticesOnCell.data(...
                1:idx1,cidxn(i)))*180/pi;

         verts(1:idx1)=ncvert.verticesOnCell.data(...
                                     1:idx1,cidxn(i));
         cells(1:idx1)=cidxn(i);

         lat(idx1+1)=NaN;
         lon(idx1+1)=NaN;

         verts(idx1+1)=NaN;
         cells(idx1+1)=cidxn(i);

         if idx3==1

          lat(idx2+1:maxidx+1)=ncvert.clatitude.data(...
                ncvert.verticesOnCell.data(...
                idx2:maxidx,cidxn(i)))*180/pi;
          lon(idx2+1:maxidx+1)=ncvert.clongitude.data(...
                ncvert.verticesOnCell.data(...
                idx2:maxidx,cidxn(i)))*180/pi;

          verts(idx2+1:maxidx+1)=ncvert.verticesOnCell.data(...
                                          idx2:maxidx,cidxn(i));
          cells(idx2+1:maxidx+1)=cidxn(i);

         elseif idx3>idx2

          lat(idx2+1:idx3+1)=ncvert.clatitude.data(...
                ncvert.verticesOnCell.data(...
                idx2:idx3,cidxn(i)))*180/pi;
          lon(idx2+1:idx3+1)=ncvert.clongitude.data(...
                ncvert.verticesOnCell.data(...
                idx2:idx3,cidxn(i)))*180/pi;

          verts(idx2+1:idx3+1)=ncvert.verticesOnCell.data(...
                                          idx2:idx3,cidxn(i));
          cells(idx2+1:idx3+1)=cidxn(i);

          lat(idx3+2)=NaN;
          lon(idx3+2)=NaN;

          verts(idx3+2)=NaN;
          cells(idx3+2)=cidxn(i);

          lat(idx4+2:maxidx+2)=ncvert.clatitude.data(...
                ncvert.verticesOnCell.data(...
                idx4:maxidx,cidxn(i)))*180/pi;
          lon(idx4+2:maxidx+2)=ncvert.clongitude.data(...
                ncvert.verticesOnCell.data(...
                idx4:maxidx,cidxn(i)))*180/pi;

          verts(idx4+2:maxidx+2)=ncvert.verticesOnCell.data(...
                                   idx4:maxidx,cidxn(i));
          cells(idx4+2:maxidx+2)=cidxn(i);

          lon(maxidx+3)=lon(1);
          lat(maxidx+3)=lat(1);

          verts(maxidx+3)=verts(1);
          cells(maxidx+3)=cidxn(i);

         end

         % append to string of points
         px1=[px1; NaN; lat(:)];
         py1=[py1; NaN; lon(:)];
         pverts1=[pverts1; NaN; verts(:)];
         pcells1=[pcells1; NaN; cells(:)];

         clear lat lon verts cells idx1 idx2 idx3 idx4

        else

         idx1=find(ncvert.verticesOnCell.data(1:maxidx,cidxn(i))==...
                  cflagcoast1(k))*xflagcoast(k);

         idx2=find(ncvert.verticesOnCell.data(1:maxidx,cidxn(i))==...
                  cflagcoast2(k))*xflagcoast(k);

         xxp=sort([idx1 idx2]);

         if ~isnan(idx1) & ~isnan(idx2)

          idx1=xxp(1);
          idx2=xxp(2);

          if idx2==idx1+1

           lat(1:idx1)=ncvert.clatitude.data(...
                ncvert.verticesOnCell.data(...
                1:idx1,cidxn(i)))*180/pi;
           lon(1:idx1)=ncvert.clongitude.data(...
                ncvert.verticesOnCell.data(...
                1:idx1,cidxn(i)))*180/pi;
           
           verts(1:idx1)=ncvert.verticesOnCell.data(1:idx1,cidxn(i));
           cells(1:idx1)=cidxn(i);

           lat(idx1+1)=NaN;
           lon(idx1+1)=NaN;

           verts(idx1+1)=NaN;
           cells(idx1+1)=cidxn(i);

           lat(idx2+1:maxidx+1)=ncvert.clatitude.data(...
                ncvert.verticesOnCell.data(...
                idx2:maxidx,cidxn(i)))*180/pi;
           lon(idx2+1:maxidx+1)=ncvert.clongitude.data(...
                ncvert.verticesOnCell.data(...
                idx2:maxidx,cidxn(i)))*180/pi;
           
           verts(idx2+1:maxidx+1)=ncvert.verticesOnCell.data(...
                                          idx2:maxidx,cidxn(i));
           cells(idx2+1:maxidx+1)=cidxn(i);

           lon(maxidx+2)=lon(1);
           lat(maxidx+2)=lat(1);

           verts(maxidx+2)=verts(1);
           cells(maxidx+2)=cidxn(i);

          else

           lat(1:maxidx)=ncvert.clatitude.data(...
                ncvert.verticesOnCell.data(...
                1:maxidx,cidxn(i)))*180/pi;
           lon(1:maxidx)=ncvert.clongitude.data(...
                ncvert.verticesOnCell.data(...
                1:maxidx,cidxn(i)))*180/pi;

           verts(1:maxidx)=ncvert.verticesOnCell.data(...
                                       1:maxidx,cidxn(i));
           cells(1:maxidx)=cidxn(i);

          end

          % append to string of points
          px2=[px2; NaN; lat(:)];
          py2=[py2; NaN; lon(:)];
          pverts2=[pverts2; NaN; verts(:)];
          pcells2=[pcells2; NaN; cells(:)];

          clear lat lon verts cells idx1 idx2

         end

        end

       else

        lat(1:maxidx)=ncvert.clatitude.data(...
                ncvert.verticesOnCell.data(...
                1:maxidx,cidxn(i)))*180/pi;
        lon(1:maxidx)=ncvert.clongitude.data(...
                ncvert.verticesOnCell.data(...
                1:maxidx,cidxn(i)))*180/pi;

        verts(1:maxidx)=ncvert.verticesOnCell.data(...
                                     1:maxidx,cidxn(i));
        cells(1:maxidx)=cidxn(i);

        lat(maxidx+1)=lat(1);
        lon(maxidx+1)=lon(1);

        verts(maxidx+1)=verts(1);
        cells(maxidx+1)=cidxn(i);

        % append to string of points
        px3=[px3; NaN; lat(:)];
        py3=[py3; NaN; lon(:)];
        pverts3=[pverts3; NaN; verts(:)];
        pcells3=[pcells3; NaN; cells(:)];

        clear lat lon verts cells

       end

      end

      % compile list of all lines
      tpx=[px0; NaN; px1; NaN; px2; NaN; px3];
      tpy=[py0; NaN; py1; NaN; py2; NaN; py3];
      tpverts=[pverts0; NaN; pverts1; NaN; pverts2; NaN; pverts3];
      tpcells=[pcells0; NaN; pcells1; NaN; pcells2; NaN; pcells3];

else

 error('no coast found')

end

% now pass sequence of vertices to sperate only 
% by one NaN
tpverts(isnan(tpx))=NaN;
tpcells(isnan(tpx))=NaN;

vlats=[];
vlons=[];
verts=[];
cells=[];
k=1;
if ~isnan(tpx(1))
 vlats(k)=tpx(1);
 vlons(k)=tpy(1);
 verts(k)=tpverts(1);
 cells(k)=tpcells(1);
 k=k+1;
end
for i=2:length(tpx)-1
 if ~(isnan(tpx(i-1)) & isnan(tpx(i)))  & ...
    ~(isnan(tpx(i-1)) & ~isnan(tpx(i)) & isnan(tpx(i+1)))
  vlats(k)=tpx(i);
  vlons(k)=tpy(i);
  verts(k)=tpverts(i);
  cells(k)=tpcells(i);
  k=k+1;
 end
end
if ~isnan(tpx(end))
 vlats(k)=tpx(end);
 vlons(k)=tpy(end);
 verts(k)=tpverts(end);
 cells(k)=tpcells(end);
end

tpx=vlats;
tpy=vlons;
tpverts=verts;
tpcells=cells;

clear vlats vlons verts cells;

k=1;
for i=2:length(tpx)
 if ~(isnan(tpx(i-1)) & isnan(tpx(i))) 
  vlats(k)=tpx(i);
  vlons(k)=tpy(i);
  verts(k)=tpverts(i);
  cells(k)=tpcells(i);
  k=k+1;
 end
end

% boundary condition on previous loop
if isnan(vlats(end-1)) & ~isnan(vlats(end))
 vlats=vlats(1:end-2);
 vlons=vlons(1:end-2);
 verts=verts(1:end-2);
 cells=cells(1:end-2);
end

tpx=[NaN vlats NaN];
tpy=[NaN vlons NaN];
tpverts=[NaN verts NaN];
tpcells=[NaN cells NaN];

% create netcdf structure

nc.attributes.title='Threshold Edge Definition';

nc.npoints.data=[1:length(tpx)];
nc.npoints.long_name='number of points on threshold';
nc.npoints.dimension={'npoints'};

nc.latitude.data=tpx;
nc.latitude.long_name='latitude of threshold';
nc.latitude.dimension={'npoints'};

nc.longitude.data=tpy;
nc.longitude.long_name='longitude of threshold';
nc.longitude.dimension={'npoints'};

nc.cells.data=tpcells;
nc.cells.long_name='cell indices';
nc.cells.dimension={'npoints'};

nc.vertices.data=tpverts;
nc.vertices.long_name='vertex indices';
nc.vertices.dimension={'npoints'};

% plot data if no output argument is specified
if nargout==0

 [x,y,z,phi,theta]=ridgepack_satfwd(nc.latitude.data,...
                                    nc.longitude.data,...
                                    centlat,centlon,...
                                    2*horizon,1.0001);
 plot3(x,y,z,'Color','m',...
       'LineWidth',0.5-sin(deg2rad(horizon))*0.4)

 hold on

 axis off

 axis equal
 view([0 0 0.4])
 axis tight

else

 nc=ridgepack_struct(nc);

end

