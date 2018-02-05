% read ns3 file and denoise tDCS artifacts
clear all; close all; clc
addpath(genpath('../../COMMON'));
nsxFilename='../data/BRA_AL/NSP2-129-256channels with micro/zHMmGZ_20171208-101841-001.ns5'; % cognitive
%nevFilename='../data/BRA_AL/NSP2-129-256channels with micro/zHMmGZ_20171208-101841-001.nev'; % cognitive

%%
fs=30000;
cellIndx=2;
channelsKeep=1:8;
%channelsKeep=1:16; % C and TM (F are in white matter)
flHz=300;
fhHz=1000;
dsr=round(fs/fhHz*0.5); % oversample by factor of 2 for resolving spikes
fsr=fs/dsr;
butterOrder=4;
dimFilter=2;

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
% GSVD parameters
cleanBlockIndx=3; % right before stim
dirtyBlockInds=4:6;

%%
% electrode subsets
cInds=1:8;
bInds=9:16;
tmInds=17:24;

%%
% block conditions
preBlockInds=1:3;
stimBlockInds=4:6;
postBlockInds=7:9;

%%
% firing rate kernel
sigma=0.1;
t=linspace(-0.5,0.5,fs);
kernel=exp(-1/(2*pi)*(t/sigma).^2);

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
data=(downsample(data.',dsr)).';
%%
% downsample?


%%
% epoch data into all blocks
nBlocks=numel(blockStartTimesStr);
blockData=zeros(size(data,1),round(blockDurationSamples/dsr),nBlocks);
for b=1:nBlocks
    blockData(:,:,b)=epochData(data,recordStartTimeStr,blockStartTimesStr{b},blockDurationSecs,fsr);
end

%%
% % denoise using GSVD
% nDirtyBlocks=numel(dirtyBlockInds);
% for b=1:nDirtyBlocks
%     [U,V,XX,C,S] = gsvd(blockData(:,:,cleanBlockIndx).',blockData(:,:,dirtyBlockInds(b)).',0);
%     gsvs=diag(C)./diag(S);
%     figure; plot(gsvs,'*k'); % use this to determine K
%     title(b);
%     K(b)=input('Enter the value of K');
%     %TODO: make a hard threshold
%     blockData(:,:,dirtyBlockInds(b))=(V*S(:,K(b)+1:end)*XX(:,K(b)+1:end).').';
% end
% %%


% figure;
% hold on
% plot(dataStim(1,:));
% plot(dataStimCleaned(1,:));

%%
% spike detection
threshStds=5;
mua=zeros(size(blockData));
for b=1:nBlocks
  mua(:,:,b) = detectSpikes(blockData(:,:,b),threshStds);
end

%%
% % convert to firing rate with gaussian kernel
% 
% clear data
% clear blockData
% 
% firing=filter(kernel,1,mua,[],2);

%%
% separate into C and TM electrodes
% 
% dsr=10;
% firingDown=downsample(firing,dsr,2);
% 
% cFiring=firing(cInds,:,:);
% tmFiring=firing(tmInds,:,:);
% 
% cFiringPre=mean(cFiring(:,:,preInds),3);
% cFiringStim=mean(cFiring(:,:,stimInds),3);
% cFiringPost=mean(cFiring(:,:,postInds),3);



% [muaStimCleaned,spikeIndsStim] = detectSpikes(dataStimCleaned,threshStds);
% [muaPost,spikeIndsPost] = detectSpikes(dataPost,threshStds);

%%
% count spikes
nSpikes=squeeze(sum(mua,2));
spikeRate=nSpikes/blockDurationSecs;

% nSpikes_c=squeeze(sum(mua(cInds,:,:),2));
% spikeRate_c=nSpikes_c/blockDurationSecs;
% 
% nSpikes_b=squeeze(sum(mua(bInds,:,:),2));
% spikeRate_b=nSpikes_b/blockDurationSecs;
% 
% nSpikes_tm=squeeze(sum(mua(tmInds,:,:),2));
% spikeRate_tm=nSpikes_tm/blockDurationSecs;


%%
% figure out tdcs onset for plot
% [~, ~, ~, H0, MN0, S0] = datevec(tdcsStartTimeStr);
% tdcsStartTimeSecs=H0*3600+MN0*60+S0;
% 
% [~, ~, ~, H1, MN1, S1] = datevec(tdcsStopTimeStr);
% tdcsStopTimeSecs=H1*3600+MN1*60+S1;


%%

% figure
% subplot(311); hold on
% plot(1:9,mean(spikeRate,1),'--ko','MarkerFaceColor','k');
% set(gca,'XTickLabel',{});
% title('C electrodes','FontWeight','normal')


% yl=ylim;
% area(linspace(3.3,6.3,2),[yl(2) yl(2)],'FaceAlpha',0.2,'BaseValue',yl(1));
% ylim(yl);
% subplot(312); hold on
% plot(1:9,mean(spikeRate_b,1),'--ko','MarkerFaceColor','k');
% set(gca,'XTickLabel',{});
% yl=ylim;
% area(linspace(3.3,6.3,2),[yl(2) yl(2)],'FaceAlpha',0.2,'BaseValue',yl(1));
% ylim(yl);
% ylabel('Mean Firing Rate (Hz)');
% title('TM electrodes','FontWeight','normal')
% subplot(313); hold on
% plot(1:9,mean(spikeRate_tm,1),'--ko','MarkerFaceColor','k');
% yl=ylim;
% area(linspace(3.3,6.3,2),[yl(2) yl(2)],'FaceAlpha',0.2,'BaseValue',yl(1));
% ylim(yl);
% xlabel('Time (hh:mm:ss)');
% title('F electrodes','FontWeight','normal')
% set(gca,'XTickLabel',blockStartTimesStr);

% figure
% subplot(311); hold on
% plot(1:9,mean(spikeRate_c,1),'--ko','MarkerFaceColor','k');
% set(gca,'XTickLabel',{});
% title('C electrodes','FontWeight','normal')
% yl=ylim;
% area(linspace(3.3,6.3,2),[yl(2) yl(2)],'FaceAlpha',0.2,'BaseValue',yl(1));
% ylim(yl);
% subplot(312); hold on
% plot(1:9,mean(spikeRate_b,1),'--ko','MarkerFaceColor','k');
% set(gca,'XTickLabel',{});
% yl=ylim;
% area(linspace(3.3,6.3,2),[yl(2) yl(2)],'FaceAlpha',0.2,'BaseValue',yl(1));
% ylim(yl);
% ylabel('Mean Firing Rate (Hz)');
% title('TM electrodes','FontWeight','normal')
% subplot(313); hold on
% plot(1:9,mean(spikeRate_tm,1),'--ko','MarkerFaceColor','k');
% yl=ylim;
% area(linspace(3.3,6.3,2),[yl(2) yl(2)],'FaceAlpha',0.2,'BaseValue',yl(1));
% ylim(yl);
% xlabel('Time (hh:mm:ss)');
% title('F electrodes','FontWeight','normal')
% set(gca,'XTickLabel',blockStartTimesStr);
% 
% print -dpng ../figures/BRA_AL_cognitive_mean_firing_rate
%%

% cInds=1:8;
% bInds=9:16;
% tmInds=17:24;
% figure;
% subplot(311);
% bar([nSpikesPre(cInds) nSpikesStimCleaned(cInds) nSpikesPost(cInds)]);
% subplot(312);
% bar([nSpikesPre(bInds) nSpikesStimCleaned(bInds) nSpikesPost(bInds)]);
% subplot(313);
% bar([nSpikesPre(tmInds) nSpikesStimCleaned(tmInds) nSpikesPost(tmInds)]);

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


