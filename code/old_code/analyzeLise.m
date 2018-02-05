clear all; close all; clc
% things to do:
% - check for bad channels (e.g. no 6 Hz power)

write=0;
outFilename='../data/Extracted tDCS/allE 6Hz.mat';
figFilename='../figures/allE 6Hz';

% for each intracerebral stimulation site, get the name of the two adjacent electrodes that were stimulated
leadFieldChannelsFilename='../data/Extracted/leadFieldChannels.mat';

% the actual data file recorded during DC stimulation
trcFilenames={'../data/Extracted tDCS/tDCS-SIN-G1-2-trial 2.trc'; ...
    '../data/Extracted tDCS/tDCS-SIN_G9-10 - trial 2.TRC'; ...
    '../data/Extracted tDCS/tDCS-SIN_G16-17.trc'; ...
    '../data/Extracted tDCS/tDCS-SIN-P1-2.TRC'; ...
    '../data/Extracted tDCS/tDCS-SIN-P8-9.TRC'};
indxTargets=[21;22;24;7;9];
strTargets={'G1-G2','G9-G10','G16-G17','P1-P2','P8-P9'};


if numel(trcFilenames)~=numel(indxTargets)
    error('length of trcFilenames must match length of indxTargets');
end
nTargets=numel(indxTargets);

fs=2048;
fHigh=5;
show=0;
nLeadFieldChannels=52;
delx=3.5;  % spacing between adjacent contacts in mm

allE=zeros(nLeadFieldChannels,nTargets);

for t=1 %:nTargets
    
    [brainDataWindow,channelNames]=readStimulationData(trcFilenames{t},fHigh,show);

    durSamples=size(brainDataWindow,2);
    %nfft=2^nextpow2(durSamples);
    nfft=durSamples;
    freqs=(0:nfft-1)/nfft*fs;
    [~,indx6Hz]=min(abs(freqs-6));
    fftBrainDataWindow=fft(brainDataWindow,nfft,2);
    %Pss=abs(fftBrainDataWindow);
    %Phi_ss=angle(fftBrainDataWindow);
    
    P2=abs(fftBrainDataWindow/durSamples);
    P1=P2(:,1:durSamples/2+1);
    P1(:,2:end-1)=2*P1(:,2:end-1);
    Pss=P1;
    %Phi_ss=angle(P1);
    

    
    %%
    % this block for trying to consider the effect of phase (probably not what
    % we want)
    % move back to time domain but only at 6 Hz
    % fft_only6=zeros(size(fftBrainDataWindow));
    % fft_only6(:,indx6Hz)=fftBrainDataWindow(:,indx6Hz);
    % brainDataWindow_only6=real(ifft(fft_only6,nfft,2));
    
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
        
        % compute electric field using magnitude only
        E(e)=(Pss(indxChannel1,indx6Hz)-Pss(indxChannel2,indx6Hz))/(delx);
        
        % if you want to consider phase
        %E(e)=mean(brainDataWindow_only6(indxChannel1)-brainDataWindow_only6(indxChannel2));
    end
    
    allE(:,t)=E;
    %%
    figure(1);
    subplot(2,3,t);
    stem(E); hold on
    stem(indxTargets(t),E(indxTargets(t)),'r');
    title(strTargets{t});
    
    
end

%%
if write
    save(outFilename,'allE','strTargets','indxTargets');
    print('-dpng',figFilename)
end

%%
% figure;
% stem(Phi_ss(:,indx6Hz));

% figure
% stem(db(Pss(:,indx6Hz)));

%%
chPwr=mean(brainDataWindow.^2,2);
figure
stem(db(chPwr));