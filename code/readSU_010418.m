% read ns3 file and denoise tDCS artifacts
clear all; close all; clc
addpath(genpath('../../COMMON'));
nsxFilename='../data/BRA_AL/NSP2-129-256channels with micro/zHMmGZ_20171208-101841-001.ns5'; % cognitive
%nevFilename='../data/BRA_AL/NSP2-129-256channels with micro/zHMmGZ_20171208-101841-001.nev'; % cognitive

%%
cellIndx=2;
channelsKeep=1:24;
flHz=300;
fhHz=1000;
butterOrder=4;
dimFilter=2;

%%
% from Laurent's email
BLOCK_DURATION_SECS=70;
TIME_START='10:23:08.000'; 
TIME_BLOCK1='10:23:16'; TIME_BLOCK2='10:24:42'; TIME_BLOCK3='10:26:06'; 
TIME_BLOCK4='10:30:24'; TIME_BLOCK5='10:32:13'; TIME_BLOCK6='10:33:43'; 
TIME_BLOCK7='10:36:16'; TIME_BLOCK8='10:37:57'; TIME_BLOCK9='10:39:18'; 

%start familiar face block1 : 10h23min16s
%end familiar face block1 : 10h24min26s
%start familiar face block2 : 10h24min42s
%end familiar face block2 : 10h25min52s
%start familiar face block3 : 10h26min06s
%end familiar face block3: 10h27min16s

%tDCS ON : 10h27min45s

%start familiar face block4 : 10h30min24s
%end familiar face block4 : 10h32min02s
%start familiar face block5: 10h32min13s
%end familiar face block5 : 10h33min23s
%start familiar face block6 : 10h33min43s
%end familiar face block6: 10h34min53s

%tDCS OFF: 10h35min38s 

%start familiar face block7: 10h36min16s
%end familiar face block7: 10h37min27s
%start familiar face block8: 10h37min57s
%end familiar face block8: 10h39min07s
%start familiar face block9: 10h39min18s
%end familiar face block9: 10h40min28s

%% 
% read in data and get array size and sampling rate
NSx = openNSx(nsxFilename);
fs=NSx.MetaTags.SamplingFreq;
data=NSx.Data{cellIndx}(channelsKeep,:);
durationMinutes=size(data,2)/fs/60
%channelLabels=[NSx.ElectrodesInfo.Label];

%%
% convert data to uV
maxAnalogs=[NSx.ElectrodesInfo(channelsKeep).MaxAnalogValue];
minAnalogs=[NSx.ElectrodesInfo(channelsKeep).MinAnalogValue];
maxDigitals=[NSx.ElectrodesInfo(channelsKeep).MaxDigiValue];
minDigitals=[NSx.ElectrodesInfo(channelsKeep).MinDigiValue];

maxAnalogs=double(maxAnalogs);
minAnalogs=double(minAnalogs);
maxDigitals=double(maxDigitals);
minDigitals=double(minDigitals);

conversion= (maxAnalogs-minAnalogs)./(maxDigitals-minDigitals);

data=double(data).*repmat(conversion.',1,size(data,2));

%%
% prefiltering
fl=flHz/(fs/2);
fh=fhHz/(fs/2);
[b,a]=butter(butterOrder,[fl fh]);
data=filter(b,a,data,[],dimFilter);


%%
% epoch data into 3 cognitive blocks
blockDurationSamples=round(BLOCK_DURATION_SECS*fs);

%dataPre = epochData(data,TIME_START,TIME_BLOCK1,BLOCK_DURATION_SECS,fs);
%dataStim = epochData(data,TIME_START,TIME_BLOCK4,BLOCK_DURATION_SECS,fs);
%dataPost = epochData(data,TIME_START,TIME_BLOCK7,BLOCK_DURATION_SECS,fs);

dataPre = epochData(data,TIME_START,TIME_BLOCK2,BLOCK_DURATION_SECS,fs);
dataStim = epochData(data,TIME_START,TIME_BLOCK5,BLOCK_DURATION_SECS,fs);
dataPost = epochData(data,TIME_START,TIME_BLOCK8,BLOCK_DURATION_SECS,fs);


%%
% denoise using GSVD
[U,V,XX,C,S] = gsvd(dataPre.',dataStim.',0);
gsvs=diag(C)./diag(S);
%figure; plot(gsvs,'*k'); % use this to determine K
%%
K=1;
dataStimCleaned=(V*S(:,K+1:end)*XX(:,K+1:end).').';

figure;
hold on
plot(dataStim(1,:));
plot(dataStimCleaned(1,:));

%%
% spike detection
threshStds=5;
[muaPre,spikeIndsPre] = detectSpikes(dataPre,threshStds);
[muaStimCleaned,spikeIndsStim] = detectSpikes(dataStimCleaned,threshStds);
[muaPost,spikeIndsPost] = detectSpikes(dataPost,threshStds);

%%
% count spikes
nSpikesPre=sum(muaPre,2);
nSpikesStimCleaned=sum(muaStimCleaned,2);
nSpikesPost=sum(muaPost,2);

cInds=1:8;
bInds=9:16;
tmInds=17:24;
figure;
subplot(311);
bar([nSpikesPre(cInds) nSpikesStimCleaned(cInds) nSpikesPost(cInds)]);
subplot(312);
bar([nSpikesPre(bInds) nSpikesStimCleaned(bInds) nSpikesPost(bInds)]);
subplot(313);
bar([nSpikesPre(tmInds) nSpikesStimCleaned(tmInds) nSpikesPost(tmInds)]);

%%
% figure; hold on
% plot(dataDetect(1,:));
% stem(mua(1,:).*dataDetect(1,:));
% plot(thresh(1)*ones(1,size(dataDetect,2)),'--');
% xlim([1 1*60*fs]);










% this section for visualizing data
% %%
% % look at impulse response of band-pass filter
% % ir=filter(b,a,[1 zeros(1,299)]);
% % figure;
% % plot(ir);
% %%
% %
% channelIndx=1;
% timeIndx=10;
% windowLength=5*fs;
% window=hamming(windowLength);
% [s,f,t]=spectrogram(data(channelIndx,:),window,[],[],fs);
% 
% %%
% % plot spectrogram
% figure;
% subplot(211);
% plot(f,abs(s(:,timeIndx)));
% subplot(212);
% imagesc([t(1) t(end)]/60,[f(1) f(end)],abs(s));
% ylim([0 2000]);
% %imagesc(abs(s));
% 
% %%
% % plot time series of one channel
% channelIndx=1;
% figure;
% plot(data(channelIndx,:));


