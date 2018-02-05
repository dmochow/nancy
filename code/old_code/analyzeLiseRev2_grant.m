clear all; close all; clc

% only consider the stimulated intracerebral pairs
% render in 3D

dataFromDisk=1;
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
precomputedFilenames={'../precomputed/tDCS-SIN-G1-2-trial 2.mat'; ...
    '../precomputed/tDCS-SIN_G9-10 - trial 2.mat'; ...
    '../precomputed/tDCS-SIN_G16-17.mat'; ...
    '../precomputed/tDCS-SIN-P1-2.mat'; ...
    '../precomputed/tDCS-SIN-P8-9.mat'};


electrodeCoordinateFilename='../data/Extracted/electrodeCoordinates.mat';
load(electrodeCoordinateFilename,'locs','labels'); % locs and labels now in memory
nLocs=size(locs,1);

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

for t=1:nTargets
    
    %[brainDataWindow,channelNames]=readStimulationData(trcFilenames{t},fHigh,show);
    trcFilename=trcFilenames{t};
    indTarget=indxTargets(t);
    precomputedFilename=precomputedFilenames{t};
    
    
    if ~dataFromDisk % select the range yourself
        [brainDataWindow,channelNames]=readStimulationData(trcFilename,fHigh,show);
        save(precomputedFilename,'brainDataWindow','channelNames');
        % save the data here so that you don't have to keep selecting the range
    else % grab epoched data from file
        load(precomputedFilename,'brainDataWindow','channelNames')
    end
    
    
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
        E(e)=0.001*(Pss(indxChannel1,indx6Hz)-Pss(indxChannel2,indx6Hz))/(delx);
        
        % try to locate this pair in 3D
        channel1strNoSpaces=strrep(channel1str,' ','');
        tmp1 = strcmp(labels, channel1strNoSpaces);
        locindx1=find(tmp1);
        channel2strNoSpaces=strrep(channel2str,' ','');
        tmp2 = strcmp(labels, channel2strNoSpaces);
        locindx2=find(tmp2);
        
        if ~isempty(locindx1) && isempty(locindx2)
            plotlocs(e,:)=locs(locindx1,:);
        elseif isempty(locindx1) && ~isempty(locindx2)
            plotlocs(e,:)=locs(locindx2,:);
        elseif ~isempty(locindx1) && ~isempty(locindx2)
            plotlocs(e,:)=mean([locs(locindx1,:);locs(locindx2,:)],1);
        else
            plotlocs(e,:)=[NaN NaN NaN];
        end
    end
    
    allE(:,t)=E;
    %%
    
    if E(indxTargets(t))<0
        E=-E;
    end
    load ../data/redgreyblue
    
    
    figure(1);
    subplot(2,3,t);
    stem(E); hold on
    stem(indxTargets(t),E(indxTargets(t)),'r');
    title(strTargets{t});
    
    figure(2);
    subplot(2,3,t);
    S=40;
    scatter3(plotlocs(:,1),plotlocs(:,2),plotlocs(:,3),S,E,'filled'); hold on
    scatter3(plotlocs(indxTargets(t),1),plotlocs(indxTargets(t),2),plotlocs(indxTargets(t),3),2*S,E(indxTargets(t)),'filled','d');
    title(strTargets{t});
    caxis([min(E) max(E)]);
    colormap(cm);
    colorbar
    axis equal
    axis off
    
end

%%
if write
    save(outFilename,'allE','strTargets','indxTargets');
    print('-dpng',figFilename)
end

