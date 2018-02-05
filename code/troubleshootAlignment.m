%%
% attempt to align NSx and TRC files
clear all; close all; clc
addpath(genpath('../../COMMON'));
nsxFilename='../data/BRA_AL/NSP2-129-256channels with micro/zHMmGZ_20171208-101841-001.ns5'; % cognitive
trcFilename='../data/BRA_AL/micromed/EEG_2014.TRC';

%%
fs=30000;
cellIndx=2;
channelsKeep=1:16; % C and TM (F are in white matter)

%%
% from Laurent's email
blockDurationSecs=70;
recordStartTimeStr='10:23:08.000'; 
blockStartTimesStr={'10:23:16';'10:24:42';'10:26:06';... %pre tDCS
'10:30:24'; '10:32:13'; '10:33:43'; ... % during tDCS
'10:36:16'; '10:37:57'; '10:39:18'}; % post tDCS
blockDurationSamples=round(blockDurationSecs*fs);
tdcsStartTimeStr='10:27:45';
tdcsStopTimeStr='10:35:38';


%%
%
NSx = openNSx(nsxFilename);
fs=NSx.MetaTags.SamplingFreq;
nsxData=NSx.Data{cellIndx}(channelsKeep,:);

% read micromed
[trcData,output]  = icem_read_micromed_trc(trcFilename);
trcData=cell2mat(trcData);
fstrc=output.SR;

%%
figure;
plot((0:size(trcData,2)-1)/fstrc,trcData(1,:));
%xlim([565600 583300]);
%timeTrc=ginput(1);

%%
figure
plot((0:size(nsxData,2)-1)/fs,nsxData(1,:));
%xlim([8195000 8292000]);
%timeNsx=ginput(1);

%%
% timeTrc=timeTrc(1)/fstrc;
% timeNsx=timeNsx(1)/fs;
% timeTrc-timeNsx

% %%
% secsElapsed = timeStr2Sample(tdcsStartTimeStr,recordStartTimeStr,fstrc)
% timeTrc
% 
% %%
% tdcsStartTimeStr_nsx = timeSample2Str(recordStartTimeStr,timeNsx,fs)

%%
% cross-check start times of blocks
startTimesSecs=double([output.NOTE_DATA.Starting_Sample])/fstrc;
timeStr=cell(5,1);
for i=1:5
    timeStr{i} = timeSample2Str(recordStartTimeStr,startTimesSecs(i),fstrc);
end
timeStr