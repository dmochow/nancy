clear all; close all; clc

csvFilename='../data/PAT_1/HUE_LIS.csv';
saveFilename='../data/Extracted/electrodeCoordinates.mat';

fid=fopen(csvFilename);
N=170;
C=textscan(fid, '%s%f%f%f%f%f%f%f%f%f%f%f%f%s%s%s',N, 'delimiter',',');
locs=[C{2} C{3} C{4}];
labels=C{1};
save(saveFilename,'locs','labels');

