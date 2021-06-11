function [seaiceheatNh,seaiceheatSh]=ridgepack_E3SM_sea_ice_cpl_heat_check(filename)

%% Initialize variables.
%filename = '/Users/afroberts/data/MODEL/E3SM/20210316/flux/monthly.txt';
%filename = '/Users/afroberts/data/MODEL/E3SM/20210316/flux/monthly_heat.txt';
%filename = '/Users/afroberts/data/MODEL/E3SM/20210316/oldlogs/flux/monthly.txt';

%% Format for each line of text:
%   column1: categorical (%C)
%	column2: double (%f)
%   column3: double (%f)
%	column4: double (%f)
%   column5: double (%f)
%	column6: double (%f)
%   column7: double (%f)
% For more information, see the TEXTSCAN documentation.
formatSpec = '%16C%15f%15f%15f%15f%15f%15f%15f%f%[^\n\r]';

%% Open the text file.
fileID = fopen(filename,'r');

%% Read columns of data according to the format.
% This call is based on the structure of the file used to generate this
% code. If an error occurs for a different file, try regenerating the code
% from the Import Tool.
dataArray = textscan(fileID, formatSpec, 'Delimiter', '', 'WhiteSpace', '', 'TextType', 'string',  'ReturnOnError', false);

%% Close the text file.
fclose(fileID);

%% Post processing for unimportable data.
% No unimportable data rules were applied during the import, so no post
% processing code is included. To generate code which works for
% unimportable data, select unimportable cells in a file and regenerate the
% script.

%% Create output variable
temp = table(dataArray{1:end-1}, 'VariableNames', {'descriptor','atm','lnd','rof','ocn','iceNh','iceSh','glc','sum'});

%% Clear temporary variables
clearvars filename formatSpec fileID dataArray ans;

%% PLOT DATA

time=[1:length(temp.sum)];
for i=1:length(temp.sum); 
 seaiceheatNh(i)=temp.iceNh(i);
 seaiceheatSh(i)=temp.iceSh(i);
end

