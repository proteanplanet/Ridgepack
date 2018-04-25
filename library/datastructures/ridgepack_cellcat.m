function sc=ridgepack_cellcat(S,charinsert)

% ridgepack_cellcat - Concatenates string elements in a cell array to single string
%
% function sc=ridgepack_cellcat(S,charinsert)
%
% This function creates a single string from many individual strings 
% in a cell array S, with elements blank-space separated.
%
% INPUT:
%
% S  - Cell array containing multiple strings. If S is a string only, 
%      then just the value of the string is returned.
%
% charinsert - Inserts a character between each string. For example,
%              you may wish to set charinsert="_" to place an underbar 
%              between the characters. This is optional, and if omitted, 
%              the default is to place a blanks space strings.
%
%
% OUTPUT:
%
% sc - Concatenated string from S
%
% Example:
% if S={'This','is','a','test'} then the output
% will be sc=['This is a test']
% 
% Ridgepack Version 1.0
% Andrew Roberts, Naval Postgraduate School, March 2018 (afrobert@nps.edu)
% 

global debug;
if debug; disp(['Entering ',mfilename,'...']); end

sc=[];

if nargin<2
 charinsert=' ';
elseif ~ischar(charinsert)
 error('charinsert must be a character string')
end

if iscell(S)
 for i=1:length(S)
	if ~ischar(S{i})
 	 error('Entry for cell array is not a string')
	end
        if strcmp(charinsert,' ') | i>1
	 sc=[sc,charinsert,char(S{i})];
        else
	 sc=[char(S{i})];
        end
 end
elseif ischar(S)
        if strcmp(charinsert,' ') 
	 sc=[charinsert,S];
        else
         sc=S;
        end
else
 error('Input is not a cell array nor a string')
end

if debug; disp(['...Leaving ',mfilename]); end
