%%
% (1) model refractory period
% 
% read ns3 file and denoise tDCS artifacts
clear all; close all; clc
addpath(genpath('../../COMMON'));
nsxFilename='../data/BRA_AL/NSP2-129-256channels with micro/zHMmGZ_20171208-101841-001.ns5'; % cognitive
trcFilename='../data/BRA_AL/micromed/EEG_2014.TRC';

%%
fs=30000;
cellIndx=2;
channelsKeep=1:16; % C and TM (F are in white matter)
flHz=300;
fhHz=3000;
butterOrder=4;
dimFilter=2;

%%
% spike detection threshold (multiple of median(abs)/0.6745)
threshStds=5;

%%
%nFlips=numel(Times);
%fps=120; % screen refresh rate
fflicker=2;
nStim=140; % number of images in block
famIndx=5:5:nStim;
nfamIndx=setdiff(1:nStim,famIndx);
timesFamSecs=(famIndx-1)/fflicker;
timesNfamSecs=(nfamIndx-1)/fflicker;
timesFamSamples=round(timesFamSecs*fs);
timesNfamSamples=round(timesNfamSecs*fs);
%%
% from Laurent's email
blockDurationSecs=70;
%recordStartTimeStr='10:23:08.000'; 
nsxDelay=3.3045; % in seconds
recordStartTimeStr='10:23:11.3045'; 
% this is the corrected time after examining the raw trc and nsx files, and
% locating the time of the giant tdcs onset artifact
%recordStartTimeStr='10:23:08.000'; % original time provided by Laurent
blockStartTimesStr={'10:23:16';'10:24:42';'10:26:06';... %pre tDCS
'10:30:24'; '10:32:13'; '10:33:43'; ... % during tDCS
'10:36:16'; '10:37:57'; '10:39:18'}; % post tDCS
% blockStartTimesStr={'10:23:14.943';'10:24:42';'10:26:06';... %pre tDCS
% '10:30:51.19'; '10:32:13'; '10:33:43'; ... % during tDCS
% '10:36:17.19'; '10:37:57'; '10:39:18'}; % post tDCS
blockDurationSamples=round(blockDurationSecs*fs);
tdcsStartTimeStr='10:27:45';
tdcsStopTimeStr='10:35:38';

% from trc file
%timeStr =
%    {'10:23:14.943' }
%    {'10:27:44.973' }
%    {'10:30:51.19.5'}
%    {'10:35:38.963' }
%    {'10:36:17.19.5'}

% PREPROCESSING %

%%
% GSVD parameters
cleanBlockIndx=3; % right before stim
dirtyBlockInds=4:6;
Kmin=0.2;

%%
% electrode subsets
cInds=1:8;
tmInds=9:16;
fInds=17:24;

%%
% block conditions
preBlockInds=1:3;
stimBlockInds=4:6;
postBlockInds=7:9;
blockGroups={preBlockInds;stimBlockInds;postBlockInds};


%%
% for epoching responses to individual images
epochDurationSecs=0.5; % to make this longer, need to handle edges of data 
epochDurationSamples=epochDurationSecs*fs;

%%
% gross artifact removal parameters
Q=50; % how many standard deviations above median(abs)/0.6745 to call it artifact
xtent=0.05*fs; % how many adjacent samples to remove alongside artifacts

%% 
% read in data and get array size and sampling rate
NSx = openNSx(nsxFilename);
fs=NSx.MetaTags.SamplingFreq;
data=NSx.Data{cellIndx}(channelsKeep,:);
durationMinutes=size(data,2)/fs/60;
nChannels=size(data,1);
channelLabels=[NSx.ElectrodesInfo.Label];

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
data=(filtfilt(b,a,data.')).';

%%
% epoch data into all blocks
nBlocks=numel(blockStartTimesStr);
blockData=zeros(size(data,1),blockDurationSamples,nBlocks);
for b=1:nBlocks
    blockData(:,:,b)=epochData(data,recordStartTimeStr,blockStartTimesStr{b},blockDurationSecs,fs);
end

%% 
% remove gross artifacts
for b=1:nBlocks
    tmpIn=blockData(:,:,b);
    tmpOut=zeroBadSamples(blockData(:,:,b),Q,xtent);
    blockData(:,:,b)=tmpOut;
end

%%
% denoise using GSVD
nDirtyBlocks=numel(dirtyBlockInds);
refData=[blockData(:,:,1) blockData(:,:,2) blockData(:,:,3)]; 
%noisyData=[blockData(:,:,4) blockData(:,:,5) blockData(:,:,6)]; 
for b=4:6
    blockData(:,:,b) = gsvdDenoise(blockData(:,:,b),refData,Kmin);
end

%%
% spike detection
mua=zeros(size(blockData));
for b=1:nBlocks
    [spikes,waveforms,mua(:,:,b),thresh] = detectSpikes(blockData(:,:,b),threshStds);
    
%     if b>3 & b<7
%         %[spikes,waveforms,mua(:,:,b),thresh] = detectSpikes(blockData(:,:,b),threshStds,fs,'neg',allThresh(:,3));
%         [spikes,waveforms,mua(:,:,b),thresh] = detectSpikes(blockData(:,:,b),threshStds,fs,'neg',mean(allThresh(:,1:3),2));
%     end
    
    allSpikes(:,b)=spikes.';
    allWaveforms(:,b)=waveforms.';
    allThresh(:,b)=thresh;
end

%%
% throw out bad spikes
for ch=1:nChannels
    for b=1:nBlocks
        [ch,b]
        spikesIn=allSpikes{ch,b};
        waveformIn=allWaveforms{ch,b};
        thresh=0.5;
        [cleanSpikes,cleanWaveform]=removeBadSpikes(spikesIn,waveformIn,thresh);
        allSpikes{ch,b}=cleanSpikes;
        allWaveforms{ch,b}=cleanWaveform;
    end
end

%%
%convert to mua
mua=zeros(nChannels,size(blockData,2),nBlocks);
for ch=1:nChannels
    for b=1:nBlocks
        mua(ch,allSpikes{ch,b},b)=1;
    end
end

% ANALYSIS %

%%
% compute PSTH for all units and blocks
allPsths=cell(nChannels,nBlocks);
binWidthMs=500;
for ch=1:nChannels
    for b=1:nBlocks
        [psth,binCenters] = computePsth(squeeze(mua(ch,:,b)).',binWidthMs,fs);
        allPsths{ch,b}=psth;
    end
end


%%
% epoch into familiar versus non-familiar
famMua=zeros(nChannels,epochDurationSamples,nBlocks,numel(timesFamSecs));
nfamMua=zeros(nChannels,epochDurationSamples,nBlocks,numel(timesNfamSecs));
for b=1:nBlocks
    for t=1:numel(timesFamSecs)
        famMua(:,:,b,t)=mua(:,timesFamSamples(t)+1:timesFamSamples(t)+epochDurationSamples,b);
        nfamMua(:,:,b,t)=mua(:,timesNfamSamples(t)+1:timesNfamSamples(t)+epochDurationSamples,b);
    end
end

%%
% create psth for each electrode and condition
binWidthMs=25;
famPsths=[];nfamPsths=[];
for ch=1:nChannels
    for b=1:nBlocks
        tmua=squeeze(famMua(ch,:,b,:));
        tmua=tmua(:,:);
        [famPsths(:,b,ch),binCenters] = computePsth(tmua,binWidthMs,fs);
        
        tmua=squeeze(nfamMua(ch,:,b,:));
        tmua=tmua(:,:);
        nfamPsths(:,b,ch) = computePsth(tmua,binWidthMs,fs);    
    end
end

%%
ylims=[0 1];
barWidth=1;
chans2plot=[3 6 7 12 15];
nChans2plot=numel(chans2plot);
blocks2plot=[1:9];
nBlocks2plot=numel(blocks2plot);
figure;
for ch=1:nChans2plot
    for b=1:nBlocks2plot
    
    hs(ch,b)=subplot(nChans2plot,nBlocks2plot,(ch-1)*nBlocks2plot+b);
    try
        bar(binCenters,famPsths(:,blocks2plot(b),chans2plot(ch)),barWidth);
    catch
        bar(squeeze(famPsths(:,blocks2plot(b),chans2plot(ch))),barWidth);
    end
    ylim(ylims);
    if ch==1, title('Familiar','FontWeight','normal'); end
   
%     hs(ch,2)=subplot(8,2,ch*2);
%     try
%         bar(binCenters,nfamPsths(:,blocks2plot,chans2plot(ch)),barWidth);
%     catch
%         bar(squeeze(nfamPsths(:,blocks2plot,chans2plot(ch))),barWidth);
%     end
%     ylim(ylims);
%     if ch==1, title('Non-Familiar','FontWeight','normal'); end
    end
end
% set(get(hs(8,1),'Xlabel'),'String','Time (s)');
% set(get(hs(1,1),'Ylabel'),'String','PSTH');
% 
% hlg=legend('Pre','tDCS','Post');
% set(hlg,'box','off');
% set(hlg,'Orientation','horizontal');
% lgPos=get(hlg,'Position');
% set(hlg,'Position',[lgPos(1) lgPos(2)-0.075 lgPos(3) lgPos(4)]);




% FIGURES %

%%
% total number of spikes in each block
grossSpikeCount=cellfun(@numel,allSpikes);
activeIndx=find(mean(grossSpikeCount,2)>300);
figure;
hs(1)=subplot(211);
%plot(grossSpikeCount(activeIndx,:).'/blockDurationSecs,'--o');
bar(grossSpikeCount(activeIndx,:).'/blockDurationSecs);
ylabel('Spike rate (Hz)');
hlg=legend({'C3','C6','C7','TM4','TM7'});
set(hlg,'orientation','horizontal');
set(hlg,'box','off');
set(hlg,'location','northwest');
ylim([0 17]);
hs(2)=subplot(212);
%plot(mean(grossSpikeCount(activeIndx,:),1)/blockDurationSecs,'--o');
bar(mean(grossSpikeCount(activeIndx,:),1)/blockDurationSecs);
ylabel('Spike rate (Hz)');

set(hs(1),'Xticklabel',{'Pre 1','Pre 2','Pre 3','TDCS 1','TDCS 2','TDCS 3','Post 1','Post 2','Post 3'});
set(hs(2),'Xticklabel',{'Pre 1','Pre 2','Pre 3','TDCS 1','TDCS 2','TDCS 3','Post 1','Post 2','Post 3'});

print -dpng ../figures/JEA_SI_spike_rate_whole_block
%%
% nRows=numel(activeIndx);
% nCols=9;
% figure;
% for r=1:nRows
%     for c=1:nCols
%         subplot(nRows,nCols,(r-1)*nCols+c);
%         bar(binCenters,allPsths{activeIndx(r),c});
%         ylim([0 20]);
%     end
% end
        

% 
% 
% %%


