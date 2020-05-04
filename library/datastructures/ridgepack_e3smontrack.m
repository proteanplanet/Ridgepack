function [cell,oncell,celllist,oncelllist,celldist,cellarea]=...
              ridgepack_e3smontrack(ncvert,tracklat,tracklon)

idx=find(ncvert.latCell.data*180/pi>60);

tic
for i=1:length(tracklat)

 [cell(i),vert,tvert,incell(i),celldist(i),vdist,cidx,vidx,tidx]=...
  ridgepack_e3smtriangulate(ncvert,tracklat(i),tracklon(i),idx);

end
toc

% distill to only cells where points fall inside
oncell=cell(incell);

% create second list removing multiple reference to same cell
k=1;
celllist(1)=cell(1);
for i=2:length(cell)
 if cell(i)~=cell(i-1)
  k=k+1;
  celllist(k)=cell(i);
 end
end

k=1;
oncelllist(1)=oncell(1);
for i=2:length(oncell)
 if oncell(i)~=oncell(i-1)
  k=k+1;
  oncelllist(k)=oncell(i);
 end
end

cellarea=ncvert.areaCell.data(cell);

