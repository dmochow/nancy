clear all; close all; clc

ROAST_PATHNAME='/Users/jacekdmochowski/PROJECTS/TACS/code/roast';
addpath(genpath(ROAST_PATHNAME));

NII_FILENAME='../data/BRA_AL/anat.nii';

%recipes={}
roast(NII_FILENAME,{'Fp1',1,'Iz',-1});
