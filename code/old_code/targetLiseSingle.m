clear all; close all; clc

%%
toSave=0;
leadFieldFilename='../../data/HUE_LI/Extracted/leadFields.mat';
locFilename='../../data/HUE_LI/HUE_LI_EEG-SEEG2.elc';
Imax=0.002;
nScalpElectrodes=24;
nBrainElectrodes=52;
load('../../data/HUE_LI/Extracted/recordingNamesExtended.mat'); % variable recordingNames now in memory
K=9;

% load colormaps
load('../../data/HUE_LI/redgreyblue.mat');
cmrgb=cm;
load('../../data/HUE_LI/redgrey.mat');
cmrg=cm;
cmusa=jmaColors('usadarkblue');


%% for this session only
targetIndx=9;
%outFilename=['../data/Extracted/montages_' num2str(K) 'electrodes_1ma/montage_' recordingNames{targetIndx} '_' num2str(K) 'electrodes.txt'];
%figFilename=['../data/Extracted/montages_' num2str(K) 'electrodes_1ma/montage_' recordingNames{targetIndx} '_' num2str(K) 'electrodes.png'];
saveFilename=['../data/../HUE_LI/Extracted/montages_' num2str(K) 'electrodes_1ma/montage_' recordingNames{targetIndx} '_' num2str(K) '_field_and_montage.mat'];


%%
load(leadFieldFilename);
V=leadFields(:,targetIndx);
head.R=leadFields;
Il1 = reciprocateL1(V,head,Imax);
El1=leadFields.'*Il1;

%%
E=zeros(52,1);
E(targetIndx)=1;
[Ilsq,optval] = solve_lsq_l1(E,head.R,Imax);
Elsq=head.R.'*Ilsq;

%%
Ilsq_down = reduceMontage(E,Il1,head.R,K,Imax);
Ilsq_down=Ilsq_down(1:end-1);
Elsq_down=head.R.'*Ilsq_down;

%%
eloc = readlocs( locFilename );
eloc_no_ref=eloc(2:end);
figure;
subplot(331)
topoplot(leadFields(:,targetIndx),eloc_no_ref,'electrodes','labels');
colormap(cmusa)
subplot(332)
topoplot_dc(Il1,eloc_no_ref,'electrodes','labels');
colormap(cmusa)
title(recordingNames{targetIndx});
subplot(333);
stem(El1); hold on
stem(targetIndx,El1(targetIndx),'r')

subplot(335)
topoplot_dc(Ilsq,eloc_no_ref,'electrodes','labels');
colormap(cmusa)
subplot(336);
stem(Elsq); hold on
stem(targetIndx,Elsq(targetIndx),'r')

subplot(338)
topoplot_dc(Ilsq_down,eloc_no_ref,'electrodes','labels');
colormap(cmusa)
subplot(339);
stem(Elsq_down); hold on
stem(targetIndx,Elsq_down(targetIndx),'r')

%% write solution to file
if toSave
    save(saveFilename,'targetIndx','El1','Elsq','Elsq_down','leadFields','Il1');
end

