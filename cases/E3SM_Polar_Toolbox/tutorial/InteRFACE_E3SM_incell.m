% InteRFACE_E3SM_incell: Tutorial Script
%
% This script demonstrates how to use the basic search function
% for determining the nearest E3SM scalar grid point, and if the
% search location is within an E3SM cell or not. 
%
% InteRFACE Tutorial, April 2020
% Andrew Roberts, LANL, afroberts@lanl.gov

clf
clear

% CHANGE THESE THINGS TO MAKE THIS WORK FOR YOU: 

% plot location
plotloc='/Users/afroberts/work/tutorial';

% grab the E3SM grid data from the gridloc
gridloc='/Users/afroberts/data/MODEL/E3SM/DECK/grid';
gridfile='E3SM_LR_V1_grid.nc';

% EVERYTHING BELOW THIS LINE SHOULD JUST WORK

cd(gridloc)

% read in grid data 
ncvert=ridgepack_clone(gridfile,{'latVertex','lonVertex',...
                                'verticesOnCell','indexToCellID',...
                                'nEdgesOnCell','edgesOnCell',...
                                'cellsOnEdge','cellsOnVertex',...
                                'edgesOnVertex'});

% splice in cell and edge latitude and longitude
nccell=ridgepack_clone(gridfile,{'latCell','lonCell'});

ncedge=ridgepack_clone(gridfile,{'latEdge','lonEdge'});

ncvert.latCell=nccell.latitude;
ncvert.lonCell=nccell.longitude;
ncvert.latEdge=ncedge.latitude;
ncvert.lonEdge=ncedge.longitude;

% add in the coastline and plot it on simply Cartesian grid
[nc,SCP,SCL]=ridgepack_e3smcoastm(ncvert);

hc=plot(nc.longitude.data,nc.latitude.data,...
        ':','Color',[0.88 0.3 0],'Linewidth',0.75);

hold on

% find search lat and lon on the mesh, finding if it is in the cell
searchlat=71.055;
searchlon=156.27;

[cell,vert,tvert,incell,cdist,vdist,cidx,vidx,tidx]=...
          ridgepack_e3smtriangulate(ncvert,searchlat,searchlon);

if incell
 disp('Inside Cell')
else
 disp('Outside Cell')
end

% plot grid cell
plot(...
 (180/pi)*ncvert.longitude.data(ncvert.verticesOnCell.data([1:ncvert.nEdgesOnCell.data(cell) 1],cell)),...
 (180/pi)*ncvert.latitude.data(ncvert.verticesOnCell.data([1:ncvert.nEdgesOnCell.data(cell) 1],cell)),...
      'b-')

hold on

% plot relevant triangulation
lats=[ncvert.latitude.data(vert)*180/pi ncvert.latitude.data(tvert)*180/pi ncvert.latCell.data(cell)*180/pi];
lons=[ncvert.longitude.data(vert)*180/pi ncvert.longitude.data(tvert)*180/pi ncvert.lonCell.data(cell)*180/pi];

plot([lons lons(1)],[lats lats(1)],'--','Color',0.5*[1 1 1])

% plot triangulation points
plot(ncvert.longitude.data(vert)*180/pi,...
     ncvert.latitude.data(vert)*180/pi,'ro')
text(ncvert.longitude.data(vert)*180/pi,...
     ncvert.latitude.data(vert)*180/pi,...
    {'','Nearest Vertex'},...
    'HorizontalAlignment','Center',...
    'VerticalAlignment','top',...
    'Color','r')

plot(ncvert.lonCell.data(cell)*180/pi,...
     ncvert.latCell.data(cell)*180/pi,'ko')
text(ncvert.lonCell.data(cell)*180/pi,...
     ncvert.latCell.data(cell)*180/pi,...
    {'Nearest Cell Center',''},...
    'HorizontalAlignment','Center',...
    'VerticalAlignment','bottom',...
    'Color','k')

plot(wrapTo360(searchlon),searchlat,'s','Color',[0 0.5 0])
text(wrapTo360(searchlon),searchlat,...
    'Inside Point $\rightarrow$   ',...
    'HorizontalAlignment','right','Color',[0 0.5 0])

% Now search for a second point
searchlat=71.00;
searchlon=156.2;

[cello,verto,tverto,incell,cdist,vdist,cidx,vidx,tidx]=...
          ridgepack_e3smtriangulate(ncvert,searchlat,searchlon);

plot(wrapTo360(searchlon),searchlat,'s','Color',[0 0.5 0])
text(wrapTo360(searchlon),searchlat,...
    'Outside Point $\rightarrow$   ',...
    'HorizontalAlignment','right','Color',[0 0.5 0])

if cello==cell & tverto==tvert
 if incell
  disp('Inside Cell')
 else
  disp('Outside Cell')
 end
else
 disp('NO MATCH')
end

% Now plot Utqiagvik (previously Barrow) for reference
searchlat=71.283;
searchlon=156.783;

plot(wrapTo360(searchlon),searchlat,'+','Color',[0 0 1])
text(wrapTo360(searchlon),searchlat,...
    'Utqiagvik ',...
    'HorizontalAlignment','right','Color',[0 0 1])

% Tidy up the plot
xlim([155.5 157])
ylim([70.9 71.5])

xlabel('Longitude ($^\circ$W)')
ylabel('Latitude ($^\circ$N)')

legend(hc,'Coast','Location','SouthEast')
legend boxoff

title('Example of Triangulation of Random Locations on Unstructured E3SM Ice-Ocean 30-60km Mesh')

% print it
cd(plotloc)
ridgepack_fprint('png',mfilename,1,2)


