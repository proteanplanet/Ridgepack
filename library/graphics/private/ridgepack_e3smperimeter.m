function [vlats,vlons,verts]=ridgepack_e3smperimeter(ncvert,cidx,infill)


if nargin<3
 infill=false;
end

if nargin<2
 error('Not enough inputs')
end

debug=false;
%debug=true;

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% First, grab all sides connecting vertices shared by
% only two cells, or not shared at all, thus defining
% the edge of an area within a cidx cell cluster.

% get size of total vertices in cells within criteria
cnlength=sum(ncvert.nEdgesOnCell.data(cidx));

% now use it to set the of vertexlist
cvertexlist=zeros([cnlength 1]);

% compile list of vertices
k=1;
for i=cidx'
 cvertexlist(k:k+ncvert.nEdgesOnCell.data(i)-1)=...
  ncvert.verticesOnCell.data(1:ncvert.nEdgesOnCell.data(i),i);
 k=k+ncvert.nEdgesOnCell.data(i);
end

% sort list of data to find number of occurrences
csortexvertex=sort(cvertexlist);

% find number of c cells per vertex
disp('Doing length search')
cdixv=zeros([length(ncvert.nVertices.data) 1]);
for i=1:length(csortexvertex)
 cdixv(csortexvertex(i))=cdixv(csortexvertex(i))+1;
end

% now calculate mask of vertices at land border noting that 
% a only a vertex bording 1 or 2 cells sits on the coast.
cpassx=NaN*zeros([length(ncvert.nVertices.data) 1]);
cpassx(cdixv==1 | cdixv==2)=1;

% now find where cells have vertex bordering 2 cells
cpassy=NaN*zeros([length(ncvert.nVertices.data) 1]);
cpassy(cdixv==2)=1;

% mask out lats and longs that are not at threshold edge
ncvert.clatitude.data=cpassx.*ncvert.latitude.data;
ncvert.clongitude.data=cpassx.*ncvert.longitude.data;

% calculate mask for cells on the threshold edge and for 
% sides with 2 adjoined vertices. Flagcell indicates
% a cell in which the cell has two adjacent 2-vertex
% locations, and flags those vertices for a later search. 
% This is to stop them being joined on the coastline.
disp('Masking to threshold edge cells only')
cpasscell=zeros([length(ncvert.nCells.data(cidx)) 1]);
k=1;
cvertexflag1=[];
cvertexflag2=[];
cflagcell=[];
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

% find shared edges on cells of ice edge
[cvflag1,I]=sort(cvertexflag1);
cvflag2=cvertexflag2(I);
cfcell=cflagcell(I);

k=1;
cflagcoast1=[];
cflagcoast2=[];
ccellflag1=[];
ccellflag2=[];
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

% only pass through threshold cells
disp('Calculating the extent line')
if length(cidxn)>0

     px0=[];
     py0=[];
     pverts0=[];

     px1=[];
     py1=[];
     pverts1=[];

     px2=[];
     py2=[];
     pverts2=[];

     px3=[];
     py3=[];
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

        lat=[];
        lon=[];
        verts=[];
        idx1=[];
        idx2=[];
        idx3=[];
        idx4=[];
        idx5=[];
        idx6=[];

        idx1=find(ncvert.verticesOnCell.data(1:maxidx,cidxn(i))==...
                 cflagcoast1(k(1)));

        idx2=find(ncvert.verticesOnCell.data(1:maxidx,cidxn(i))==...
                 cflagcoast2(k(1)));

        idx3=find(ncvert.verticesOnCell.data(1:maxidx,cidxn(i))==...
                 cflagcoast1(k(2)));

        idx4=find(ncvert.verticesOnCell.data(1:maxidx,cidxn(i))==...
                 cflagcoast2(k(2)));

        idx5=find(ncvert.verticesOnCell.data(1:maxidx,cidxn(i))==...
                 cflagcoast1(k(3)));

        idx6=find(ncvert.verticesOnCell.data(1:maxidx,cidxn(i))==...
                 cflagcoast2(k(3)));

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
        lat(idx1+1)=NaN;
        lon(idx1+1)=NaN;

        verts(idx1+1)=NaN;

        lat(idx2+1:idx3+1)=ncvert.clatitude.data(...
               ncvert.verticesOnCell.data(...
               idx2:idx3,cidxn(i)))*180/pi;
        lon(idx2+1:idx3+1)=ncvert.clongitude.data(...
               ncvert.verticesOnCell.data(...
               idx2:idx3,cidxn(i)))*180/pi;

        verts(idx2+1:idx3+1)=ncvert.verticesOnCell.data(...
                                idx2:idx3,cidxn(i));

        lat(idx3+2)=NaN;
        lon(idx3+2)=NaN;

        verts(idx3+2)=NaN;

        if idx5==1

         lat(idx4+2:maxidx+2)=ncvert.clatitude.data(...
               ncvert.verticesOnCell.data(...
               idx4:maxidx,cidxn(i)))*180/pi;
         lon(idx4+2:maxidx+2)=ncvert.clongitude.data(...
               ncvert.verticesOnCell.data(...
               idx4:maxidx,cidxn(i)))*180/pi;

         verts(idx4+2:maxidx+2)=ncvert.verticesOnCell.data(...
                                         idx4:maxidx,cidxn(i));

        elseif idx5>idx4

         lat(idx4+2:idx5+2)=ncvert.clatitude.data(...
               ncvert.verticesOnCell.data(...
               idx4:idx5,cidxn(i)))*180/pi;
         lon(idx4+2:idx5+2)=ncvert.clongitude.data(...
               ncvert.verticesOnCell.data(...
               idx4:idx5,cidxn(i)))*180/pi;

         verts(idx4+2:idx5+2)=ncvert.verticesOnCell.data(...
                                         idx4:idx5,cidxn(i));

         lat(idx5+3)=NaN;
         lon(idx5+3)=NaN;

         verts(idx5+3)=NaN;

         lat(idx6+3:maxidx+3)=ncvert.clatitude.data(...
               ncvert.verticesOnCell.data(...
               idx6:maxidx,cidxn(i)))*180/pi;
         lon(idx6+3:maxidx+3)=ncvert.clongitude.data(...
               ncvert.verticesOnCell.data(...
               idx6:maxidx,cidxn(i)))*180/pi;

         verts(idx6+3:maxidx+3)=ncvert.verticesOnCell.data(...
                                          idx6:maxidx,cidxn(i));

         lon(maxidx+4)=lon(1);
         lat(maxidx+4)=lat(1);

         verts(maxidx+4)=verts(1);

        end

        % append to string of points
        px0=[px0; NaN; lat(:)];
        py0=[py0; NaN; lon(:)];
        pverts0=[pverts0; NaN; verts(:)];


       elseif size(k,2)==2 

        lat=[];
        lon=[];
        verts=[];
        idx1=[];
        idx2=[];
        idx3=[];
        idx4=[];

        idx1=find(ncvert.verticesOnCell.data(1:maxidx,cidxn(i))==...
                 cflagcoast1(k(1)));

        idx2=find(ncvert.verticesOnCell.data(1:maxidx,cidxn(i))==...
                 cflagcoast2(k(1)));

        idx3=find(ncvert.verticesOnCell.data(1:maxidx,cidxn(i))==...
                 cflagcoast1(k(2)));

        idx4=find(ncvert.verticesOnCell.data(1:maxidx,cidxn(i))==...
                 cflagcoast2(k(2)));

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

        lat(idx1+1)=NaN;
        lon(idx1+1)=NaN;

        verts(idx1+1)=NaN;

        if idx3==1

         lat(idx2+1:maxidx+1)=ncvert.clatitude.data(...
               ncvert.verticesOnCell.data(...
               idx2:maxidx,cidxn(i)))*180/pi;
         lon(idx2+1:maxidx+1)=ncvert.clongitude.data(...
               ncvert.verticesOnCell.data(...
               idx2:maxidx,cidxn(i)))*180/pi;

         verts(idx2+1:maxidx+1)=ncvert.verticesOnCell.data(...
                                         idx2:maxidx,cidxn(i));

        elseif idx3>idx2

         lat(idx2+1:idx3+1)=ncvert.clatitude.data(...
               ncvert.verticesOnCell.data(...
               idx2:idx3,cidxn(i)))*180/pi;
         lon(idx2+1:idx3+1)=ncvert.clongitude.data(...
               ncvert.verticesOnCell.data(...
               idx2:idx3,cidxn(i)))*180/pi;

         verts(idx2+1:idx3+1)=ncvert.verticesOnCell.data(...
                                         idx2:idx3,cidxn(i));

         lat(idx3+2)=NaN;
         lon(idx3+2)=NaN;

         verts(idx3+2)=NaN;

         lat(idx4+2:maxidx+2)=ncvert.clatitude.data(...
               ncvert.verticesOnCell.data(...
               idx4:maxidx,cidxn(i)))*180/pi;
         lon(idx4+2:maxidx+2)=ncvert.clongitude.data(...
               ncvert.verticesOnCell.data(...
               idx4:maxidx,cidxn(i)))*180/pi;

         verts(idx4+2:maxidx+2)=ncvert.verticesOnCell.data(...
                                  idx4:maxidx,cidxn(i));

         lon(maxidx+3)=lon(1);
         lat(maxidx+3)=lat(1);

         verts(maxidx+3)=verts(1);

        end

        px1=[px1; NaN; lat(:)];
        py1=[py1; NaN; lon(:)];
        pverts1=[pverts1; NaN; verts(:)];


       else

        lat=[];
        lon=[];
        verts=[];
        idx1=[];
        idx2=[];

        idx1=find(ncvert.verticesOnCell.data(1:maxidx,cidxn(i))==...
                 cflagcoast1(k));

        idx2=find(ncvert.verticesOnCell.data(1:maxidx,cidxn(i))==...
                 cflagcoast2(k));

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

          lat(idx1+1)=NaN;
          lon(idx1+1)=NaN;

          verts(idx1+1)=NaN;

          lat(idx2+1:maxidx+1)=ncvert.clatitude.data(...
               ncvert.verticesOnCell.data(...
               idx2:maxidx,cidxn(i)))*180/pi;
          lon(idx2+1:maxidx+1)=ncvert.clongitude.data(...
               ncvert.verticesOnCell.data(...
               idx2:maxidx,cidxn(i)))*180/pi;
          
          verts(idx2+1:maxidx+1)=ncvert.verticesOnCell.data(...
                                         idx2:maxidx,cidxn(i));

          lon(maxidx+2)=lon(1);
          lat(maxidx+2)=lat(1);

          verts(maxidx+2)=verts(1);

         else

          lat(1:maxidx)=ncvert.clatitude.data(...
               ncvert.verticesOnCell.data(...
               1:maxidx,cidxn(i)))*180/pi;
          lon(1:maxidx)=ncvert.clongitude.data(...
               ncvert.verticesOnCell.data(...
               1:maxidx,cidxn(i)))*180/pi;

          verts(1:maxidx)=ncvert.verticesOnCell.data(...
                                      1:maxidx,cidxn(i));

         end

         % append to string of points
         px2=[px2; NaN; lat(:)];
         py2=[py2; NaN; lon(:)];
         pverts2=[pverts2; NaN; verts(:)];

        end

       end

      else

       lat=[];
       lon=[];
       verts=[];

       lat(1:maxidx)=ncvert.clatitude.data(...
               ncvert.verticesOnCell.data(...
               1:maxidx,cidxn(i)))*180/pi;
       lon(1:maxidx)=ncvert.clongitude.data(...
               ncvert.verticesOnCell.data(...
               1:maxidx,cidxn(i)))*180/pi;

       verts(1:maxidx)=ncvert.verticesOnCell.data(...
                                    1:maxidx,cidxn(i));
       lat(maxidx+1)=lat(1);
       lon(maxidx+1)=lon(1);

       verts(maxidx+1)=verts(1);

       % append to string of points
       px3=[px3; NaN; lat(:)];
       py3=[py3; NaN; lon(:)];
       pverts3=[pverts3; NaN; verts(:)];

      end

     end

     % compile list of all lines
     tpx=[px0; NaN; px1; NaN; px2; NaN; px3];
     tpy=[py0; NaN; py1; NaN; py2; NaN; py3];
     tpverts=[pverts0; NaN; pverts1; NaN; pverts2; NaN; pverts3];

else

 error('no coast found')

end

% now pass sequence of vertices to sperate by one NaN

% initialize
tpverts(isnan(tpx))=NaN;
vlats=[];
vlons=[];
verts=[];
k=0;

% initial condition
if ~isnan(tpx(1))
 k=k+1;
 vlats(k)=tpx(1);
 vlons(k)=tpy(1);
 verts(k)=tpverts(1);
end

% cycle through all cases
for i=2:length(tpx)-1
 if ~(isnan(tpx(i-1)) & isnan(tpx(i))) & ...
    ~(isnan(tpx(i-1)) & ~isnan(tpx(i)) & isnan(tpx(i+1)))
  k=k+1;
  vlats(k)=tpx(i);
  vlons(k)=tpy(i);
  verts(k)=tpverts(i);
 end
end

% boundary condition
if ~isnan(tpx(end))
 k=k+1;
 vlats(k)=tpx(end);
 vlons(k)=tpy(end);
 verts(k)=tpverts(end);
end

% finalize
tpx=vlats;
tpy=vlons;
tpverts=verts;

% Now remove cases of multiple NaNs

% initialize
vlats=[];
vlons=[];
verts=[];
k=0;

% Initial condition
if ~isnan(tpx(1))
 k=k+1;
 vlats(k)=tpx(1);
 vlons(k)=tpy(1);
 verts(k)=tpverts(1);
end

% cycle through looking for consecutive NaNs
for i=2:length(tpx)
 if ~(isnan(tpx(i-1)) & isnan(tpx(i))) 
  k=k+1;
  vlats(k)=tpx(i);
  vlons(k)=tpy(i);
  verts(k)=tpverts(i);
 end
end

% boundary condition on previous loop
if isnan(vlats(end-1)) & ~isnan(vlats(end))
 vlats=vlats(1:end-2);
 vlons=vlons(1:end-2);
 verts=verts(1:end-2);
end

% finalize lats, lons and vertices
vlats=[NaN vlats NaN];
vlons=[NaN vlons NaN];
verts=[NaN verts NaN];
npoin=length(verts);

% BUG IS SOMEWHERE IN HERE (ABOVE)


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Now order vertices into closed loops by joining together
% line segments with the same bonding vertices.

% find line segments - they are seperated by NaNs
idx=find(isnan(verts));

% account for case where first is a NaN
if idx(1)==1
 idx=idx(2:end);
end

% break into line segments and number them
segment=zeros(size(verts));
segment(1:idx(1))=1;
k=2;
for i=[idx(1)+1:length(verts)];
 segment(i)=k;
 if ~isempty(find(idx(:)==i))
  k=k+1;
 end
end

% set up the flag for a segment already found
% where the first element is a NaN
segmentfound=zeros(size(verts));
segmentfound(1)=1;

% set up the start index for a segment
startidx=2;
startvertex=verts(startidx);

% begin the count of the number of close contours
linecount=1;

% now hunt through the indices of line segments, and join 
% together segments in a continuous close-loop line 
for i=idx

 if segmentfound(startidx)==0

  if debug; disp('Starting new closed loop'); end

  endidx=i-1;
  endvertex=verts(endidx);

  jdx=find(verts==endvertex);
  jdx=jdx(jdx~=endidx & segmentfound(jdx)~=1);

  kdx=find(verts==startvertex);
  kdx=kdx(kdx~=startidx & segmentfound(jdx)~=1);

  segmentfound(startidx:endidx+1)=1;

  vert=verts(startidx:i-1);
  vidx=[startidx:i-1];

  if length(kdx)>1
   error('Found line-start intersection that should not exist')
  end

  while ~isempty(jdx) & ~isempty(kdx)

   % check for intersections, first looking for 'false'
   % vertices where only a single vertex is seperated by a 
   % a NaN and then exit gracefully if such a case exists. 

   if length(jdx)>1
    if all(verts(jdx)==vert(end))
     for isdf=1:length(jdx)
      if isnan(verts(jdx(isdf)-1)) & ...
         isnan(verts(jdx(isdf)+1))
       error('loan vertex found that does not define an edge')
      end
     end
    end
    error('Found line intersection that should not exist')
   end

   % case of sequential-order vertices
   if isnan(verts(jdx+1));

    if debug; disp('sequential'); end
    %disp('sequential')

    startx=find(isnan(verts(1:jdx)));
    vert=[vert verts(jdx:-1:startx(end)+1)];
    vidx=[vidx jdx:-1:(startx(end)+1)];

    segmentfound(startx(end)+1:jdx+1)=1;
    endidx=startx(end)+1;
    endvertex=vert(end);

   % case of reversed-order vertices
   elseif isnan(verts(jdx-1));

    if debug; disp('reversed'); end
    %disp('reversed');

    endx=find(isnan(verts(jdx:end)));
    vert=[vert verts(jdx+1:jdx+endx(1)-2)];
    vidx=[vidx jdx+1:jdx+endx(1)-2];

    segmentfound(jdx:jdx+endx(1)-1)=1;
    endidx=jdx+endx(1)-2;
    endvertex=vert(end);

   % have found the middle of a segment, which 
   % theoretically should not exist
   else

    error('Found intersection of lines')

   end

   oldjdx=jdx;
   alljdx=find(verts==endvertex);
   jdx=alljdx(alljdx~=endidx & segmentfound(alljdx)~=1);

  end

  if isempty(jdx) & ~isempty(kdx) & ~isempty(find(kdx==alljdx))

   if i==idx(1)
    allvert=vert;
    alllats=vlats(vidx);
    alllons=vlons(vidx);
   else
    allvert=[allvert NaN vert];
    alllats=[alllats NaN vlats(vidx)];
    alllons=[alllons NaN vlons(vidx)];
   end

   finalvert{linecount}.vert=vert;
   finalvert{linecount}.lats=vlats(vidx);
   finalvert{linecount}.lons=vlons(vidx);

   linecount=linecount+1;

   lastvert=vert;
   lastvidx=vidx;
   lastlats=vlats(vidx);
   lastlons=vlons(vidx);

   vert=[];
   vidx=[];

  elseif isempty(jdx) & ~isempty(kdx)

   disp('ERROR-----------------')

   disp(['startvertex: ',num2str(startvertex)])
   disp(['endvertex: ',num2str(endvertex)])
   disp(['startidx: ',num2str(startidx)])
   disp(['endidx: ',num2str(endidx)])

   centlat=mean(vlats(vidx));
   centlon=mean(vlons(vidx));
   horizon=5;

   ridgepack_satview(centlat,centlon,horizon,1,2); 
 
   [x,y,z,phi,theta]=ridgepack_satfwd(vlats(vidx),...
                                       vlons(vidx),...
                                       centlat,centlon,...
                                       2*horizon,1.0001);
   plot3(x,y,z,'Color','k',...
          'LineWidth',0.5-sin(deg2rad(horizon))*0.4)

   [x,y,z,phi,theta]=ridgepack_satfwd(sdlat,...
                                       sdlon,...
                                       centlat,centlon,...
                                       2*horizon,1.0001);
   plot3(x,y,z,'Color','g',...
          'LineWidth',0.5-sin(deg2rad(horizon))*0.4)

   [x,y,z,phi,theta]=ridgepack_satfwd(sdlat2,...
                                       sdlon2,...
                                       centlat,centlon,...
                                       2*horizon,1.0001);
   plot3(x,y,z,'Color','b',...
          'LineWidth',0.5-sin(deg2rad(horizon))*0.4)

   [x,y,z,phi,theta]=ridgepack_satfwd(vlats(vidx(end)),...
                                       vlons(vidx(end)),...
                                       centlat,centlon,...
                                       2*horizon,1.0001);
 
   plot3(x,y,z,'o','Color','r')

   [x,y,z,phi,theta]=ridgepack_satfwd(ttpx,...
                                       ttpy,...
                                       centlat,centlon,...
                                       2*horizon,1.0001);
   plot3(x,y,z,'Color','b',...
          'LineWidth',0.5-sin(deg2rad(horizon))*0.4)

   disp(['vidx(end): ',num2str(verts(vidx(end)))])

   axis off

   axis equal
   view([0 0 0.4])
   axis tight

   error('The contour is not closed')

  elseif isempty(jdx) & isempty(kdx)

   error('Cannot find matching start or end index')

  end

 else

  if isempty(find(startvertex==allvert(:)))
   error('Starting vertex not found in current lines')
  else
   if debug; disp('Bypass'); end
  end

 end

 if all(segmentfound(:)==1)

  disp(['Found ',num2str(linecount-1),' closed contours'])
  break

 else

  startidx=i+1;
  startvertex=verts(startidx);

 end

end

if infill

 vlats=alllats(end:-1:1);
 vlons=alllons(end:-1:1);
 verts=allvert(end:-1:1);

else

 vlats=alllats(1:1:end);
 vlons=alllons(1:1:end);
 verts=allvert(1:1:end);

end

