%%
% (1) check that nanBadSamples is not cutting out spikes
% (2) model refractory period
% 
% read ns3 file and denoise tDCS artifacts
clear all; close all; clc
addpath(genpath('../../COMMON'));
nsxFilename='../data/BRA_AL/NSP2-129-256channels with micro/zHMmGZ_20171208-101841-001.ns5'; % cognitive
%nevFilename='../data/BRA_AL/NSP2-129-256channels with micro/zHMmGZ_20171208-101841-001.nev'; % cognitive
timingFilenames={...
    '../data/BRA_AL/timing/bra_al1/bra_al_all_sinstim_Oddball_Fam_2Hz_Fam3_2017-12-08-10h15m44s.mat',...
    '../data/BRA_AL/timing/bra_al1/bra_al_all_sinstim_Oddball_Fam_2Hz_Fam6_2017-12-08-10h17m10s.mat',...
    '../data/BRA_AL/timing/bra_al1/bra_al_all_sinstim_Oddball_Fam_2Hz_Fam1_2017-12-08-10h18m34s.mat',...
    '../data/BRA_AL/timing/bra_al2/bra_al2_all_sinstim_Oddball_Fam_2Hz_Fam6_2017-12-08-10h23m20s.mat',...
    '../data/BRA_AL/timing/bra_al2/bra_al2_all_sinstim_Oddball_Fam_2Hz_Fam3_2017-12-08-10h24m41s.mat',...
    '../data/BRA_AL/timing/bra_al2/bra_al2_all_sinstim_Oddball_Fam_2Hz_Fam1_2017-12-08-10h26m11s.mat',...
    '../data/BRA_AL/timing/bra_al3/bral_al3_all_sinstim_Oddball_Fam_2Hz_Fam1_2017-12-08-10h28m45s.mat',...
    '../data/BRA_AL/timing/bra_al3/bral_al3_all_sinstim_Oddball_Fam_2Hz_Fam6_2017-12-08-10h30m25s.mat',...
    '../data/BRA_AL/timing/bra_al3/bral_al3_all_sinstim_Oddball_Fam_2Hz_Fam3_2017-12-08-10h31m46s.mat'};

trcFilename='../data/BRA_AL/micromed/EEG_2014.TRC';

precomputedFilename='../precomputed/BRA_AL/preprocessedBlackrock/blockData.mat';

%%
fs=30000;
cellIndx=2;
%channelsKeep=1:8;
channelsKeep=1:16; % C and TM (F are in white matter)
flHz=300;
fhHz=3000;
dsr=1;
fsr=fs/dsr;
butterOrder=4;
dimFilter=2;

%%
% spike detection threshold (multiple of median(abs)/0.6745)
threshStds=4;

%%
%nFlips=numel(Times);
%fps=120; % screen refresh rate
fflicker=2;
nStim=140; % number of images in block
famIndx=5:5:nStim;
nfamIndx=setdiff(1:nStim,famIndx);
timesFamSecs=(famIndx-1)/fflicker;
timesNfamSecs=(nfamIndx-1)/fflicker;
timesFamSamples=round(timesFamSecs*fsr);
timesNfamSamples=round(timesNfamSecs*fsr);
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
blockDurationSamples=round(blockDurationSecs*fs);

tdcsStartTimeStr='10:27:45';
tdcsStopTimeStr='10:35:38';

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
epochDurationSamples=epochDurationSecs*fsr;

%%
% gross artifact removal parameters
Q=15; % how many standard deviations above mean to call it artifact
xtent=0.05*fs; % how many adjacent samples to remove alongside artifacts

%%
% read micromed
%[trcData,output]  = icem_read_micromed_trc(trcFilename);


%% 
% read in data and get array size and sampling rate
NSx = openNSx(nsxFilename);
fs=NSx.MetaTags.SamplingFreq;
data=NSx.Data{cellIndx}(channelsKeep,:);
durationMinutes=size(data,2)/fs/60
nChannels=size(data,1);
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
% debug time string from laurent (3.3 s delay between trc and nsx)
%tmpEpochedData = epochData(data,recordStartTimeStr,'10:27:45',5,fs);



%%
% prefiltering
fl=flHz/(fs/2);
fh=fhHz/(fs/2);
[b,a]=butter(butterOrder,[fl fh]);
data=(filtfilt(b,a,data.')).';

%%
% epoch data into all blocks
nBlocks=numel(blockStartTimesStr);
blockData=zeros(size(data,1),round(blockDurationSamples/dsr),nBlocks);
for b=1:nBlocks
    blockData(:,:,b)=epochData(data,recordStartTimeStr,blockStartTimesStr{b},blockDurationSecs,fsr);
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
%refData=blockData(:,:,3);
refData=[blockData(:,:,1) blockData(:,:,2) blockData(:,:,3)]; % check that length mismatch doesn't cause problems
for b=4:6
    blockData(:,:,b) = gsvdDenoise(blockData(:,:,b),refData,Kmin);
end

%%
tBlockData=blockData(:,:,2);
nfft=2^nextpow2(size(tBlockData,2));
B=fft(tBlockData,nfft,2);
Pbb=B.*conj(B);
freqs=(0:nfft-1)/nfft*fs;
figure;
plot(freqs,mean(Pbb,1)); xlim([0 100]);

%%
% spike detection
mua=zeros(size(blockData));
for b=1:nBlocks
    [spikes,waveforms,mua(:,:,b),thresh] = detectSpikes(blockData(:,:,b),threshStds);
    
    if b>3 & b<7
        %[spikes,waveforms,mua(:,:,b),thresh] = detectSpikes(blockData(:,:,b),threshStds,fs,'neg',allThresh(:,3));
        [spikes,waveforms,mua(:,:,b),thresh] = detectSpikes(blockData(:,:,b),threshStds,fs,'neg',mean(allThresh(:,1:3),2));
    end
    
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

%%
% epoch into familiar versus non-familiar
famMua=zeros(nChannels,epochDurationSamples,nBlocks,numel(timesFamSecs));
nfamMua=zeros(nChannels,epochDurationSamples,nBlocks,numel(timesNfamSecs));
for b=1:nBlocks
    for t=1:numel(timesFamSecs)
        famMua(:,:,b,t)=mua(:,timesFamSamples(t)+1:timesFamSamples(t)+epochDurationSamples,b);
        nfamMua(:,:,b,t)=mua(:,timesNfamSamples(t)+1:timesNfamSamples(t)+epochDurationSamples,b);
        
        % lfp to check timing and spike detection
        famLfp(:,:,b,t)=blockData(:,timesFamSamples(t)+1:timesFamSamples(t)+epochDurationSamples,b);
        nfamLfp(:,:,b,t)=blockData(:,timesNfamSamples(t)+1:timesNfamSamples(t)+epochDurationSamples,b);
        
    end
end

%%
% create psth for each electrode and condition
binWidthMs=50;
famPsths=[];nfamPsths=[];
for ch=1:nChannels
    for bg=1:numel(blockGroups)
        tmua=squeeze(famMua(ch,:,blockGroups{bg},:));
        tmua=tmua(:,:);
        [famPsths(:,bg,ch),binCenters] = computePsth(tmua,binWidthMs,fsr);
        
        tmua=squeeze(nfamMua(ch,:,blockGroups{bg},:));
        tmua=tmua(:,:);
        nfamPsths(:,bg,ch) = computePsth(tmua,binWidthMs,fsr);    
    end
end

%%
ylims=[0 1];
barWidth=1;
chans2plot=1:8;
nChans2plot=numel(chans2plot);
figure;
for ch=1:nChans2plot
    
    hs(ch,1)=subplot(8,2,ch*2-1);
    try
        bar(binCenters,famPsths(:,:,chans2plot(ch)),barWidth);
    catch
        bar(squeeze(famPsths(:,:,chans2plot(ch))),barWidth);
    end
    ylim(ylims);
    if ch==1, title('Familiar','FontWeight','normal'); end
   
    hs(ch,2)=subplot(8,2,ch*2);
    try
        bar(binCenters,nfamPsths(:,:,chans2plot(ch)),barWidth);
    catch
        bar(squeeze(nfamPsths(:,:,chans2plot(ch))),barWidth);
    end
    ylim(ylims);
    if ch==1, title('Non-Familiar','FontWeight','normal'); end
    
end
set(get(hs(8,1),'Xlabel'),'String','Time (s)');
set(get(hs(1,1),'Ylabel'),'String','PSTH');

hlg=legend('Pre','tDCS','Post');
set(hlg,'box','off');
set(hlg,'Orientation','horizontal');
lgPos=get(hlg,'Position');
set(hlg,'Position',[lgPos(1) lgPos(2)-0.075 lgPos(3) lgPos(4)]);


%%
% total number of spikes in each block
grossSpikeCount=cellfun(@numel,allSpikes);
figure;
subplot(121);
plot(grossSpikeCount.','--o');
subplot(122);
plot(mean(grossSpikeCount,1),'--o');


%%
tmpEpochedData = epochData(data,recordStartTimeStr,tdcsStartTimeStr,5,fs);


