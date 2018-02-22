function [outset]=ridgepack_setdiff(firstset,secondset)

% ridgepack_setdiff - Generates a cell array as the difference between two cell arrays
%
% function [outset]=ridgepack_setdiff(firstset,secondset)
%
% This function generates the difference between two sets presented as 1D
% cell arrays.  This is similar to the matlab function setdiff, except the 
% output is not sorted and occurs in the same order as elements in 
% firstset
%
% Input:
% firstset   - cell array with elements that you want to test for differences
%              from the second set
% secondeset - comparing cell array
%
% Output:
% outset - cell array containing elements that are different in firstset
%          from those in second set
%
% Example:
% firstset={'x','y','z','a'}
% secondset={'x','z'}
% outset={'y','a'}
% 
% In this example, the standard matlab setdiff would provide
% outset={'a','y'}
%
% Andrew Roberts, Naval Postgraduate School, March 2018  (afrobert@nps.edu)
%

global debug;
if debug; disp(['Entering ',mfilename,'...']); end

if ~iscell(firstset)
	error('firstset is not a cell array');
elseif ~iscell(secondset)
	error('secondset is not a cell array');
end

m=0;
outset={};
for j=1:length(firstset)
        if ~any(strcmp(char(firstset{j}),secondset))
	           m=m+1;
	           outset{m}=firstset{j};
	end
end

if debug; disp(['...Leaving ',mfilename]); end

