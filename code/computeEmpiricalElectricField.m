clear all; close all; clc
addpath(genpath('../../COMMON'));

dataFromDisk=1;
trcFilename='../data/BRA_AL/micromed/EEG_2012.TRC';
xlsFilename='../data/BRA_AL/shit from Laurent/BRA_AL_Electrodes.xls';
precomputedFilename='../precomputed/BRA_AL/intracerebral_voltages_during_6Hz_scalp_stim';
fhigh=5; % high-pass cutoff frequency
show=1; % whether to show the results of epoching the data
delx=3.5;  % spacing between adjacent contacts in mm

%%
% read in the trc file and compute the intracerebral voltages during scalp
% stimulation
if ~dataFromDisk
    [brainDataWindow,channelNames,fs]=readStimulationData(trcFilename,fhigh,show);
    save(precomputedFilename,'brainDataWindow','channelNames','fs');
else
    load(precomputedFilename,'brainDataWindow','channelNames','fs');
end

channelNamesNoSpace=channelNames;
channelNamesNoSpace=strrep(channelNamesNoSpace,' ','');


%%
% now that we have the intracerebral voltages, compute the electric field
% across the brain
nfft=2^nextpow2(size(brainDataWindow,2));
freqs=(0:nfft-1)/nfft*fs;
ftarget=6;
[~,fo]=min(abs(freqs-ftarget));
B=fft(brainDataWindow,nfft,2);

B2=abs(B/nfft); % double sided
B1=B2(:,1:nfft/2+1); % single sided
B1(:,2:end-1)=2*B1(:,2:end-1); % add negative frequencies
V=B1(:,fo);

%%
% compute the electric field here (for real)
nChannels=numel(channelNames);
refChannelIndx=zeros(nChannels,1);
channelTable=[];
for ch=1:nChannels
    [thisLetterStr,thisChannelNum] = extractChannelLetter(channelNames{ch});
    % now try to find the reference for this electrode, if it exists
    refLetterStr=thisLetterStr; refChannelNum=thisChannelNum+1;
    refChannelIndx(ch) = getChannelIndx(refLetterStr,refChannelNum,channelNames);
    
    if ~isnan(refChannelIndx(ch))
        channelTable=cat(1,channelTable,[ch refChannelIndx(ch)]);
    end
end

nChannelsField=size(channelTable,1);

% finally
E=0.001*(V(channelTable(:,1))-V(channelTable(:,2)))/(delx);  % microvolts to millivolts, then divided by mm to get mV/mm = V/m
Em=abs(E); 
%%
% read the electrode locations from spreadsheet and render electric field
[num,txt,raw]=xlsread(xlsFilename);
xlsChannelNames=txt(:,1);
xlsChannelNames=strrep(xlsChannelNames,' ','');

for ch=1:nChannelsField
    thisChannelName1=channelNamesNoSpace{channelTable(ch,1)};
    ind1=find(strcmp(xlsChannelNames,thisChannelName1),1) ;
    if ~isempty(ind1)
        try
            locs(ch,:)=num(ind1,2:4);
        catch
            locs(ch,:)=[NaN,NaN,NaN];
        end
    else
        locs(ch,:)=[NaN,NaN,NaN];
    end  
end

%%
figure
stem(Em);
set(gca,'XTick',1:nChannelsField);
set(gca,'XTickLabel',channelNames(channelTable(:,1)));
set(gca,'XTickLabelRotation',90);
xt = get(gca, 'XTick');
set(gca, 'FontSize', 6)
ylabel('Electric field magnitude (V/m)');

save('../precomputed/empiricalElectricField','Em','E','channelNames','channelTable');

%%

indMissing=isnan(mean(locs,2));
Em_missingLocs=Em;
Em_missingLocs(indMissing)=0;

figure
S=60;
scatter3(locs(:,1),locs(:,2),locs(:,3),S,Em_missingLocs,'filled'); hold on
cm=jmaColors('usa');
colormap(cm);
%scatter3(plotlocs(indxTargets(t),1),plotlocs(indxTargets(t),2),plotlocs(indxTargets(t),3),2*S,E(indxTargets(t)),'filled','d');
%title(strTargets{t});
%caxis([min(E) max(E)]);
%colormap(cm);
%colorbar
%axis equal
%axis off

%%
% for correlating with the changes during cognitive task, reshape E into
% full channel set (153)

Em_full=NaN(nChannels,1);
Em_full(channelTable(:,1))=Em;

E_full=NaN(nChannels,1);
E_full(channelTable(:,1))=E;
save('../precomputed/BRA_AL_Em_full','Em_full','E_full');