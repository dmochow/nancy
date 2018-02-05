clear all; close all; clc
addpath(genpath('../../COMMON'));
%blkFilename='../data/BRA_AL/NSP2-129-256channels with
%micro/zHMmGZ_20171208-100203-001.ns5'; % targeting
%blkFilename='../data/BRA_AL/NSP2-129-256channels with micro/zHMmGZ_20171208-101841-001.ns5'; % cognitive
nevFilename='../data/BRA_AL/NSP2-129-256channels with micro/zHMmGZ_20171208-101841-001.nev'; % cognitive
%%
%NSx = openNSx(blkFilename);
%data=NSx;
%data.Data=data.Data{1};
%%
%channels=1:24;
%filter=1500;
%Spikes = findSpikes(data,'channels',channels);

%%
output= openNEV(nevFilename);
spikeTimesSamples=output.Data.Spikes.TimeStamp;
spikeElectrodes=output.Data.Spikes.Electrode;

nElectrodes=double(max(spikeElectrodes));
nSamples=double(max(spikeTimesSamples));
mua=zeros(nElectrodes,nSamples);

spikeInds = sub2ind(size(mua),spikeElectrodes,spikeTimesSamples);
mua(spikeInds)=1;

%%
stimSpikeInds= spikeTimesSamples>6*60*30000 & spikeTimesSamples<12*60*30000;
waves=output.Data.Spikes.Waveform(:,stimSpikeInds);
for i=1:size(waves,2)
    plot(waves(:,i));
    title(i);
    pause
end



