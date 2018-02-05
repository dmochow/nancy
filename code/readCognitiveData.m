clear all; close all; clc
addpath(genpath('../../COMMON'));

dataFromDisk=0;
trcFilename='../data/BRA_AL/micromed/EEG_2014.TRC';
xlsFilename='../data/BRA_AL/shit from Laurent/BRA_AL_Electrodes.xls';
eFieldFilename='../precomputed/BRA_AL_Em_full';
%precomputedFilename='../precomputed/BRA_AL/intracerebral_voltages_during_cognitive_task';
fhigh=0.7; % high-pass cutoff frequency
show=1; % whether to show the results of epoching the data
%delx=3.5;  % spacing between adjacent contacts in mm
BLOCK_DURATION_SECS=round(4*60);
fTargets=[2 4 8];



%%
if ~dataFromDisk
    [brainDataWindow,channelNames,fs]=readStimulationData(trcFilename,fhigh,show);
    %save(precomputedFilename,'brainDataWindow','channelNames','fs');
else
    %load(precomputedFilename,'brainDataWindow','channelNames','fs');
end
nChannels=numel(channelNames);
BLOCK_DURATION_SAMPLES=round(BLOCK_DURATION_SECS*fs);


%%
TIME_BLOCK1='10:23:14.943';
TIME_BLOCK2='10:30:51,019';
TIME_BLOCK3='10:36:17.019'; 

%TIME_TDCS_ON='10:35:38.963';
%TIME_TDCS_OFF='10:36:17.019';

[~, ~, ~, H1, MN1, S1] = datevec(TIME_BLOCK1);
TIME_BLOCK_1_SECS=H1*3600+MN1*60+S1;

[~, ~, ~, H2, MN2, S2] = datevec(TIME_BLOCK2);
TIME_BLOCK_2_SECS=H2*3600+MN2*60+S2;

[~, ~, ~, H3, MN3, S3] = datevec(TIME_BLOCK3);
TIME_BLOCK_3_SECS=H3*3600+MN3*60+S3;

% offset everything to onset of block 1 (assumes this is recording onset!)
TIME_BLOCK_3_SECS=TIME_BLOCK_3_SECS-TIME_BLOCK_1_SECS;
TIME_BLOCK_2_SECS=TIME_BLOCK_2_SECS-TIME_BLOCK_1_SECS;
TIME_BLOCK_1_SECS=TIME_BLOCK_1_SECS-TIME_BLOCK_1_SECS;

TIME_BLOCK_1_SAMPLES=round(TIME_BLOCK_1_SECS*fs);
TIME_BLOCK_2_SAMPLES=round(TIME_BLOCK_2_SECS*fs);
TIME_BLOCK_3_SAMPLES=round(TIME_BLOCK_3_SECS*fs);



%%
brainDataPre=brainDataWindow(:,TIME_BLOCK_1_SAMPLES+1:TIME_BLOCK_1_SAMPLES+BLOCK_DURATION_SAMPLES);
brainDataStim=brainDataWindow(:,TIME_BLOCK_2_SAMPLES+1:TIME_BLOCK_2_SAMPLES+BLOCK_DURATION_SAMPLES);
brainDataPost=brainDataWindow(:,TIME_BLOCK_3_SAMPLES+1:TIME_BLOCK_3_SAMPLES+BLOCK_DURATION_SAMPLES);


%%
nfft=2^nextpow2(BLOCK_DURATION_SAMPLES);
freqs=(0:nfft-1)/nfft*fs;
B1=fft(brainDataPre,nfft,2);
B2=fft(brainDataStim,nfft,2);
B3=fft(brainDataPost,nfft,2);

P1=B1.*conj(B1);
P2=B2.*conj(B2);
P3=B3.*conj(B3);

%%
nTargets=numel(fTargets);
fos=zeros(nTargets,1);
for t=1:nTargets
    [~,fo]=min(abs(freqs-fTargets(t)));
    fos(t)=fo;
end

powersPre=mean(P1(:,fos),2);
powersStim=mean(P2(:,fos),2);
powersPost=mean(P3(:,fos),2);

delta1=(powersStim-powersPre)./powersPre;
delta2=(powersPost-powersPre)./powersPre;



%%
% correlate changes with the field strength
load(eFieldFilename,'Em_full','E_full');
nonNanInds=~isnan(Em_full);
[r,p]=corrcoef(delta1(nonNanInds),Em_full(nonNanInds));
[r2,p2]=corrcoef(delta1(nonNanInds),E_full(nonNanInds));
[r3,p3]=corrcoef(powersStim(nonNanInds),Em_full(nonNanInds));
[r4,p4]=corrcoef(powersStim(nonNanInds),E_full(nonNanInds));
[r5,p5]=corrcoef(powersPost(nonNanInds),Em_full(nonNanInds));
[r6,p6]=corrcoef(powersPost(nonNanInds),E_full(nonNanInds));
[r7,p7]=corrcoef(powersPre(nonNanInds),Em_full(nonNanInds));
[r8,p7]=corrcoef(powersPre(nonNanInds),E_full(nonNanInds));

%%

% show differences as a function of electrode
figure
stem(delta1);
set(gca,'XTick',1:nChannels);
set(gca,'XTickLabel',channelNames);
set(gca,'XTickLabelRotation',90);
xt = get(gca, 'XTick');
set(gca, 'FontSize', 6)

figure
stem(delta2);
set(gca,'XTick',1:nChannels);
set(gca,'XTickLabel',channelNames);
set(gca,'XTickLabelRotation',90);
xt = get(gca, 'XTick');
set(gca, 'FontSize', 6)

%%
% figure;
% plot(freqs,mean(abs(B1),1));
% xlim([0 50]);

