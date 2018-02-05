clear all; close all; clc

% things to do: 
% - account for influence of phase (flips)
% - check for bad channels (e.g. no 6 Hz power)

% for each intracerebral stimulation site, get the name of the two adjacent electrodes that were stimulated
leadFieldChannelsFilename='../data/Extracted/leadFieldChannels.mat'; 

% the actual data file recorded during DC stimulation
%trcFilename='../data/Extracted tDCS/tDCS-SIN_G9-10 - trial 2.TRC';
%trcFilename='../data/Extracted tDCS/tDCS-DC_G16-17.trc';
%trcFilename='../data/Extracted tDCS/tDCS-SIN-G1-2-trial 2.trc';
trcFilename='../data/Extracted tDCS/tDCS-DC-P1-2.TRC';
%trcFilename='../data/Extracted tDCS/tDCS-DC-P8-9.TRC';

fs=2048;
fHigh=0;
show=0;
nLeadFieldChannels=52;
delx=3.5;  % spacing between adjacent contacts in mm

%%
[brainDataWindow,channelNames]=readStimulationData(trcFilename,fHigh,show);
%brainDataWindow=brainDataWindow-repmat(mean(brainDataWindow,1),size(brainDataWindow,1),1);
%% 
durSamples=size(brainDataWindow,2);
nfft=2^nextpow2(durSamples);
freqs=(0:nfft-1)/nfft*fs;
[~,indxDC]=min(abs(freqs-0));
fftBrainDataWindow=fft(brainDataWindow,nfft,2);
%Pss=fftBrainDataWindow.*conj(fftBrainDataWindow);
Pss=abs(fftBrainDataWindow);

%% convert electric potentials to electric field strength
load(leadFieldChannelsFilename);
E=zeros(nLeadFieldChannels,1);
for e=1:nLeadFieldChannels
    channel1str=leadFieldChannel1{e};
    channel2str=leadFieldChannel2{e};
    tmp1 = strcmp(channelNames, channel1str);
    tmp2 = strcmp(channelNames, channel2str);
    indxChannel1=find(tmp1);
    indxChannel2=find(tmp2);
    E(e)=(Pss(indxChannel1,indxDC)-Pss(indxChannel2,indxDC))/(delx);
end
% %%
figure;
stem(E);
% 
% %%
% [~,sortIndx]=sort(Pss(:,indx6Hz),'descend');
% channelNames(sortIndx);