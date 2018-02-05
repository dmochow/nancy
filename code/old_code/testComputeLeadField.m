clear all; close all; clc
trcFilenames={'../data/PAT_1/EEG_1.TRC','../data/PAT_1/EEG_3.TRC','../data/PAT_1/EEG_5.TRC'...
    '../data/PAT_1/EEG_7.TRC','../data/PAT_1/EEG_9.TRC','../data/PAT_1/EEG_11.TRC'};
triggerChannels=[81,85,90,66,70,76]; % needs to be adapted per file, unfortunately
locFilename='../data/PAT_1/HUE_LI_EEG-SEEG2.elc';
trc2elc=[11 12 13 14 1 2 3 4 5 6 ...
    7 8 9 10 23 24 20 21 22 15 16 17 18 19];
threshold=2000;
fhigh=10;
show=0;

eloc = readlocs( locFilename );
eloc_no_ref=eloc(2:end);
titles={'B2-B3','B6-B7','B12-B13','A2-A3','A6-A7','A12-A13'};
for s=1:6
    muScalpReordered = computeLeadField(trcFilenames{s},triggerChannels(s),locFilename,trc2elc,threshold,fhigh,show);
    subplot(2,3,s);
    topoplot(muScalpReordered,eloc_no_ref,'electrodes','numbers');
    title(titles{s});
end


