%%
% (1) check that spike extraction not picking up tDCS artifact
% (2) determine appropriate firing rate kernel

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

%%
fs=30000;
cellIndx=2;
%channelsKeep=1:8;
channelsKeep=1:16; % C and TM (F are in white matter)
flHz=300;
fhHz=1000;
dsr=round(fs/fhHz*0.5); % oversample by factor of 2 for resolving spikes
fsr=fs/dsr;
butterOrder=4;
dimFilter=2;

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
tmInds=9:16;
fInds=17:24;

%%
% block conditions
preBlockInds=1:3;
stimBlockInds=4:6;
postBlockInds=7:9;

%%
% firing rate kernel
sigma=0.005;
t=linspace(-0.25,0.25,fsr);
kernel=exp(-1/(2*pi)*(t/sigma).^2);
figure; plot(t,kernel);
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

%%
% spike detection
threshStds=5;
mua=zeros(size(blockData));
for b=1:nBlocks
  mua(:,:,b) = detectSpikes(blockData(:,:,b),threshStds);
end

%%
% convert to firing rate with gaussian kernel
%clear data blockData
firing=filter(kernel,1,mua,[],2);


%%
% epoch into familiar versus non-familiar
epochDurationSecs=0.5; % to make this longer, need to handle edges of data 
epochDurationSamples=epochDurationSecs*fsr;
famFiringEpoched=zeros(size(firing,1),epochDurationSamples,nBlocks,numel(timesFamSecs));
nfamFiringEpoched=zeros(size(firing,1),epochDurationSamples,nBlocks,numel(timesNfamSecs));
for b=1:nBlocks
    for t=1:numel(timesFamSecs)
        famFiringEpoched(:,:,b,t)=firing(:,timesFamSamples(t)+1:timesFamSamples(t)+epochDurationSamples,b);
        nfamFiringEpoched(:,:,b,t)=firing(:,timesNfamSamples(t)+1:timesNfamSamples(t)+epochDurationSamples,b);
    end
end

muFamFiringEpoched=mean(famFiringEpoched,4);
muNfamFiringEpoched=mean(nfamFiringEpoched,4);

pwrFamFiringEpoched=muFamFiringEpoched.^2;
pwrNfamFiringEpoched=muNfamFiringEpoched.^2;

tmp1=mean(pwrFamFiringEpoched(:,:,preBlockInds),3);
tmp2=mean(pwrFamFiringEpoched(:,:,stimBlockInds),3);
tmp3=mean(pwrFamFiringEpoched(:,:,postBlockInds),3);

tmp4=mean(pwrNfamFiringEpoched(:,:,preBlockInds),3);
tmp5=mean(pwrNfamFiringEpoched(:,:,stimBlockInds),3);
tmp6=mean(pwrNfamFiringEpoched(:,:,postBlockInds),3);

%%
figure
for b=1:9
    subplot(3,3,b);
    plot(muFamFiringEpoched(:,:,b).');
end

figure
for b=1:9
    subplot(3,3,b);
    plot(muNfamFiringEpoched(:,:,b).');
end


%%




% %%
% muFiringPre=mean(firing(:,:,preBlockInds),3); 
% muFiringStim=mean(firing(:,:,stimBlockInds),3); 
% muFiringPost=mean(firing(:,:,postBlockInds),3); 
% figure; 
% subplot(311); 
% imagesc(muFiringPre);
% subplot(312);
% imagesc(muFiringStim);
% subplot(313);
% imagesc(muFiringPost);





