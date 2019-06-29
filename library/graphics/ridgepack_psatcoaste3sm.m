function [nc]=ridgepack_psatcoaste3sm(ncvert,cgrid,coastname,...
                                      centlat,centlon,horizon)

if nargin<2
 cgrid==true;
elseif ~islogical(cgrid)
 error('cgrid must be either true or false')
end

if nargin<3
 coastname='E3SM_ice_ocean_coast';
elseif ~ischar(coastname)
 error('local must be either true or false')
end

if nargin<4
 centlat=90;
elseif nargin>3 & nargout==1
 error('you are asking to both plot and produce coastal data')
elseif ~isnumeric(centlat)
 error('centlat should be a real number between -90 and 90')
elseif centlat>90 | centlat<-90
 error('centlat should be between -90 and 90')
end

if nargin<5
 centlon=0;
elseif ~isnumeric(centlon)
 error('centlon should be a real number between -90 and 90')
elseif centlat>180 | centlat<-180
 error('centlon should be between -180 and 180')
end

if nargin<6
 horizon=90;
elseif ~isnumeric(horizon)
 error('horizon must be a real number between 0 and 90')
elseif horizon<0 | horizon>90
 error('horizon must be between 0 and 90')
end

% get size of total vertices in cells
nlength=sum(ncvert.nEdgesOnCell.data(:));

% now use it to set the of vertexlist
vertexlist=zeros([nlength 1]);

% compile list of vertices
k=1;
for i=1:length(ncvert.nCells.data)
 vertexlist(k:k+ncvert.nEdgesOnCell.data(i)-1)=...
  ncvert.verticesOnCell.data(1:ncvert.nEdgesOnCell.data(i),i);
 k=k+ncvert.nEdgesOnCell.data(i);
end

% sort list of data to find number of occurrences
sortexvertex=sort(vertexlist);

% find number of cells per vertex
disp('Doing length search')
dixv=zeros([length(ncvert.nVertices.data) 1]);
for i=1:length(sortexvertex)
 dixv(sortexvertex(i))=dixv(sortexvertex(i))+1;
end

% now calculate mask of vertices at land border 
% noting that a only a vertex bording 1 or 2 cells
% sits on the coast.
passx=NaN*zeros([length(ncvert.nVertices.data) 1]);
passx(dixv==1 | dixv==2)=1;

% now find where cells have vertex bordering 2 cells
passy=NaN*zeros([length(ncvert.nVertices.data) 1]);
passy(dixv==2)=1;

% mask out lats and longs that are not on a coast
ncvert.latitude.data=passx.*ncvert.latitude.data;
ncvert.longitude.data=passx.*ncvert.longitude.data;

% mask out areas not in the coastal area
if nargout==0
 disp('Reducing to just the plotting area')
 maxth=deg2rad(horizon);
 for i=1:length(ncvert.nCells.data)

  maxidx=ncvert.nEdgesOnCell.data(i);

  la=ncvert.latitude.data(ncvert.verticesOnCell.data(1:maxidx,i));
  lo=ncvert.longitude.data(ncvert.verticesOnCell.data(1:maxidx,i));

  [x,y,z,ph,th]=ridgepack_satfwd(rad2deg(la),rad2deg(lo),...
                                centlat,centlon,horizon,1);

  % filter cells no in frame, and find cropping limit
  if all(isnan(x))
   ncvert.latitude.data(ncvert.verticesOnCell.data(1:maxidx,i))=NaN;
   ncvert.longitude.data(ncvert.verticesOnCell.data(1:maxidx,i))=NaN;
  end

 end
end

% calculate mask for cells on the coast and for 
% sides with 2 adjoined vertices. Flagcell indicates
% a cell in which the cell has two adjacent 2-vertex
% locations, and flags those vertices for a later
% search. This is to stop them being joined on the 
% coastline.
disp('Masking to coastal cells only')
passcell=zeros([length(ncvert.nCells.data) 1]);
k=1;
for i=1:length(ncvert.nCells.data)

 if any(passx(ncvert.verticesOnCell.data...
       (1:ncvert.nEdgesOnCell.data(i),i))==1)
  passcell(i)=1;
 end

 for j=1:ncvert.nEdgesOnCell.data(i)

  if j<ncvert.nEdgesOnCell.data(i) & ...
     all(passy(ncvert.verticesOnCell.data(j:j+1,i))==1)
   flagcell(k)=i;
   vflag=sort(ncvert.verticesOnCell.data(j:j+1,i));
   vertexflag1(k)=vflag(1);
   vertexflag2(k)=vflag(2);
   k=k+1;
  elseif j==ncvert.nEdgesOnCell.data(i) & ...
     all(passy(ncvert.verticesOnCell.data([1 j],i))==1)
   flagcell(k)=i;
   vflag=sort(ncvert.verticesOnCell.data([1 j],i));
   vertexflag1(k)=vflag(1);
   vertexflag2(k)=vflag(2);
   k=k+1;
  end

 end

end
idxn=find(passcell==1);

% find shared edges on coastal cells
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

% only pass through coastal cells
disp('Calculating the coastal line')
if length(idxn)>0

      px0=[];
      py0=[];

      px1=[];
      py1=[];

      px2=[];
      py2=[];

      px3=[];
      py3=[];

      for i=1:length(idxn)

       maxidx=ncvert.nEdgesOnCell.data(idxn(i));

       if (any(cellflag1==idxn(i)) | any(cellflag2==idxn(i))) & cgrid

        if any(cellflag1==idxn(i))
         k=find(cellflag1==idxn(i));
        else
         k=find(cellflag2==idxn(i));
        end

        if size(k,2)==4 & cgrid

         error('Encountered four shared coastal sides')

        elseif size(k,2)==3 & cgrid

         idx1=find(ncvert.verticesOnCell.data(1:maxidx,idxn(i))==...
                  flagcoast1(k(1)));

         idx2=find(ncvert.verticesOnCell.data(1:maxidx,idxn(i))==...
                  flagcoast2(k(1)));

         idx3=find(ncvert.verticesOnCell.data(1:maxidx,idxn(i))==...
                  flagcoast1(k(2)));

         idx4=find(ncvert.verticesOnCell.data(1:maxidx,idxn(i))==...
                  flagcoast2(k(2)));

         idx5=find(ncvert.verticesOnCell.data(1:maxidx,idxn(i))==...
                  flagcoast1(k(3)));

         idx6=find(ncvert.verticesOnCell.data(1:maxidx,idxn(i))==...
                  flagcoast2(k(3)));

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

         lat(1:idx1)=ncvert.latitude.data(...
                ncvert.verticesOnCell.data(...
                1:idx1,idxn(i)))*180/pi;
         lon(1:idx1)=ncvert.longitude.data(...
                ncvert.verticesOnCell.data(...
                1:idx1,idxn(i)))*180/pi;

         lat(idx1+1)=NaN;
         lon(idx1+1)=NaN;

         lat(idx2+1:idx3+1)=ncvert.latitude.data(...
                ncvert.verticesOnCell.data(...
                idx2:idx3,idxn(i)))*180/pi;
         lon(idx2+1:idx3+1)=ncvert.longitude.data(...
                ncvert.verticesOnCell.data(...
                idx2:idx3,idxn(i)))*180/pi;
         lat(idx3+2)=NaN;
         lon(idx3+2)=NaN;

         if idx5==1

          lat(idx4+2:maxidx+2)=ncvert.latitude.data(...
                ncvert.verticesOnCell.data(...
                idx4:maxidx,idxn(i)))*180/pi;
          lon(idx4+2:maxidx+2)=ncvert.longitude.data(...
                ncvert.verticesOnCell.data(...
                idx4:maxidx,idxn(i)))*180/pi;

         elseif idx5>idx4

          lat(idx4+2:idx5+2)=ncvert.latitude.data(...
                ncvert.verticesOnCell.data(...
                idx4:idx5,idxn(i)))*180/pi;
          lon(idx4+2:idx5+2)=ncvert.longitude.data(...
                ncvert.verticesOnCell.data(...
                idx4:idx5,idxn(i)))*180/pi;

          lat(idx5+3)=NaN;
          lon(idx5+3)=NaN;

          lat(idx6+3:maxidx+3)=ncvert.latitude.data(...
                ncvert.verticesOnCell.data(...
                idx6:maxidx,idxn(i)))*180/pi;
          lon(idx6+3:maxidx+3)=ncvert.longitude.data(...
                ncvert.verticesOnCell.data(...
                idx6:maxidx,idxn(i)))*180/pi;

          lon(maxidx+4)=lon(1);
          lat(maxidx+4)=lat(1);

         end

         px0=[px0; NaN; lat(:)];
         py0=[py0; NaN; lon(:)];

         clear lat lon idx1 idx2 idx3 idx4 idx5 idx6

        elseif size(k,2)==2 & cgrid

         idx1=find(ncvert.verticesOnCell.data(1:maxidx,idxn(i))==...
                  flagcoast1(k(1)));

         idx2=find(ncvert.verticesOnCell.data(1:maxidx,idxn(i))==...
                  flagcoast2(k(1)));

         idx3=find(ncvert.verticesOnCell.data(1:maxidx,idxn(i))==...
                  flagcoast1(k(2)));

         idx4=find(ncvert.verticesOnCell.data(1:maxidx,idxn(i))==...
                  flagcoast2(k(2)));

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

         lat(1:idx1)=ncvert.latitude.data(...
                ncvert.verticesOnCell.data(...
                1:idx1,idxn(i)))*180/pi;
         lon(1:idx1)=ncvert.longitude.data(...
                ncvert.verticesOnCell.data(...
                1:idx1,idxn(i)))*180/pi;

         lat(idx1+1)=NaN;
         lon(idx1+1)=NaN;

         if idx3==1

          lat(idx2+1:maxidx+1)=ncvert.latitude.data(...
                ncvert.verticesOnCell.data(...
                idx2:maxidx,idxn(i)))*180/pi;
          lon(idx2+1:maxidx+1)=ncvert.longitude.data(...
                ncvert.verticesOnCell.data(...
                idx2:maxidx,idxn(i)))*180/pi;

         elseif idx3>idx2

          lat(idx2+1:idx3+1)=ncvert.latitude.data(...
                ncvert.verticesOnCell.data(...
                idx2:idx3,idxn(i)))*180/pi;
          lon(idx2+1:idx3+1)=ncvert.longitude.data(...
                ncvert.verticesOnCell.data(...
                idx2:idx3,idxn(i)))*180/pi;

          lat(idx3+2)=NaN;
          lon(idx3+2)=NaN;

          lat(idx4+2:maxidx+2)=ncvert.latitude.data(...
                ncvert.verticesOnCell.data(...
                idx4:maxidx,idxn(i)))*180/pi;
          lon(idx4+2:maxidx+2)=ncvert.longitude.data(...
                ncvert.verticesOnCell.data(...
                idx4:maxidx,idxn(i)))*180/pi;

          lon(maxidx+3)=lon(1);
          lat(maxidx+3)=lat(1);

         end

         px1=[px1; NaN; lat(:)];
         py1=[py1; NaN; lon(:)];

         clear lat lon idx1 idx2 idx3 idx4

        else

         idx1=find(ncvert.verticesOnCell.data(1:maxidx,idxn(i))==...
                  flagcoast1(k));

         idx2=find(ncvert.verticesOnCell.data(1:maxidx,idxn(i))==...
                  flagcoast2(k));

         xxp=sort([idx1 idx2]);

         idx1=xxp(1);
         idx2=xxp(2);

         if idx2==idx1+1

          lat(1:idx1)=ncvert.latitude.data(...
                ncvert.verticesOnCell.data(...
                1:idx1,idxn(i)))*180/pi;
          lon(1:idx1)=ncvert.longitude.data(...
                ncvert.verticesOnCell.data(...
                1:idx1,idxn(i)))*180/pi;

          lat(idx1+1)=NaN;
          lon(idx1+1)=NaN;

          lat(idx2+1:maxidx+1)=ncvert.latitude.data(...
                ncvert.verticesOnCell.data(...
                idx2:maxidx,idxn(i)))*180/pi;
          lon(idx2+1:maxidx+1)=ncvert.longitude.data(...
                ncvert.verticesOnCell.data(...
                idx2:maxidx,idxn(i)))*180/pi;

          lon(maxidx+2)=lon(1);
          lat(maxidx+2)=lat(1);

         else

          lat(1:maxidx)=ncvert.latitude.data(...
                ncvert.verticesOnCell.data(...
                1:maxidx,idxn(i)))*180/pi;
          lon(1:maxidx)=ncvert.longitude.data(...
                ncvert.verticesOnCell.data(...
                1:maxidx,idxn(i)))*180/pi;

         end

         px2=[px2; NaN; lat(:)];
         py2=[py2; NaN; lon(:)];

         clear lat lon

        end

       else

        lat(1:maxidx)=ncvert.latitude.data(...
                ncvert.verticesOnCell.data(...
                1:maxidx,idxn(i)))*180/pi;
        lon(1:maxidx)=ncvert.longitude.data(...
                ncvert.verticesOnCell.data(...
                1:maxidx,idxn(i)))*180/pi;

        lon(maxidx+1)=lon(1);
        lat(maxidx+1)=lat(1);

        px3=[px3; NaN; lat(:)];
        py3=[py3; NaN; lon(:)];

        clear lat lon

       end

      end

      % compile list of all lines
      px=[px0; NaN; px1; NaN; px2; NaN; px3];
      py=[py0; NaN; py1; NaN; py2; NaN; py3];

else

 error('no coast found')

end

% create netcdf file of coast
if cgrid
 nc.attributes.title=[coastname,' E3SM Coastline on C-grid'];
else
 nc.attributes.title=[coastname,' E3SM Coastline on B-grid'];
end

nc.npoints.data=[1:length(px)];
nc.npoints.long_name='number of points';
nc.npoints.dimension={'npoints'};

nc.latitude.data=px;
nc.latitude.long_name='latitude of coast';
nc.latitude.dimension={'npoints'};

nc.longitude.data=py;
nc.longitude.long_name='longitude of coast';
nc.longitude.dimension={'npoints'};

if nargout==1
 nc=ridgepack_struct(nc);
else
 [x,y,z,phi,theta]=ridgepack_satfwd(nc.latitude.data,...
                                    nc.longitude.data,...
                                    centlat,centlon,...
                                    2*horizon,1.0001);
 plot3(x,y,z,'Color',0.25*[1 1 1],...
       'LineWidth',0.5-sin(deg2rad(horizon))*0.4)
end

