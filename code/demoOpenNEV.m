clear all; close all; clc;

% experiment 1, pre tDCS
nevFilename='../data/JEA_SI/JEA_SI_Blackrock_129-256/qiLGkaVk_20180202-102636-001.nev'; 
nsxFilename='../data/JEA_SI/JEA_SI_Blackrock_129-256/qiLGkaVk_20180202-102636-001.ns5'; 
% experiment 1, during tDCS
%nevFilename='../data/JEA_SI/JEA_SI_Blackrock_129-256/qiLGkaVk_20180202-103531-001.nev'; % 

% experiment 2, post tDCS
% ../data/JEA_SI/JEA_SI_Blackrock_129-256/qiLGkaVk_20180202-104125-001.ns5
% ../data/JEA_SI/JEA_SI_Blackrock_129-256/qiLGkaVk_20180202-104848-001.ns5
% ../data/JEA_SI/JEA_SI_Blackrock_129-256/qiLGkaVk_20180202-105922-001.ns5
% ../data/JEA_SI/JEA_SI_Blackrock_129-256/qiLGkaVk_20180202-110654-001.ns5

NEV = openNEV(nevFilename);
timeStamps=NEV.Data.SerialDigitalIO.TimeStampSec;
blockOnsetsSec=timeStamps([find(diff([0 timeStamps])>10)])
% timeStamps(1) is the time in seconds of first block start


NSx = openNSx(nsxFilename);
nsxData=NSx.Data{2};
%%
% read in micromed 
micromedFilename='../data/JEA_SI/JEA_SI_Micromed/Familiar Face Reco 2Hz before tDCS.TRC';
[data,output]  = icem_read_micromed_trc(micromedFilename);
data=data{1};
size(data,2)/output.SR

%%
figure; 
subplot(211);
plot(data(1,:));
subplot(212);
plot(nsxData(1,:));

