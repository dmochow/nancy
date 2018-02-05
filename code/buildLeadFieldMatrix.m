clear all; close all; clc
addpath(genpath('/Users/jacekdmochowski/PROJECTS/COMMON'));
%
roastDataPath='../data/BRA_AL/';
filenames=dir([roastDataPath '*result.mat']);
load(fullfile(roastDataPath,filenames(1).name));

%%
mean(isnan(ef_mag(:)))