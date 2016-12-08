function [sampleLocationsOut] = sampleLocations
% Usage:  [sampleLocationsOut] = sampleLocations;

% Required input:
% None 

% Outputs:
% A table documenting the sample locations as part of HMS

% Load the datasets
outfilename = websave ('sampleLocations.csv','https://raw.githubusercontent.com/CatchmentSci/Teaching-2016-17/master/RiverNames.csv');
websave('readtext.m', 'https://raw.githubusercontent.com/CatchmentSci/Glaisdale-Beck-diversion-scheme/master/readtext.m'); % Download dependancy
[sampleLocationsOut,~] = readtext(outfilename, ',', '','','textual'); % read in the comma delimeted data

