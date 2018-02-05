clear all; close all; clc
dataPath='../data/extracted/';
trcFilenames=dir([dataPath '*.trc']);
nRecords=numel(trcFilenames);
nScalpElectrodes=24;
%%
load('../data/Extracted/triggerStringsExtended.mat'); % variable triggerStrings now in memory
load('../data/Extracted/recordingNamesExtended.mat'); % variable recordingNames now in memory
%%
locFilename='../data/PAT_1/HUE_LI_EEG-SEEG2.elc';
trc2elc=[11 12 13 14 1 2 3 4 5 6 ...
    7 8 9 10 23 24 20 21 22 15 16 17 18 19];
threshold=2000;
fhigh=10;
show=0;

%% loop through all recordings
badElecIndx=[];
leadFields=zeros(nScalpElectrodes,nRecords);
for r=1:nRecords
    try
        thisTrcFilename=fullfile(dataPath,trcFilenames(r).name);
        
        % find trigger channel
        [~,output]  = icem_read_micromed_trc(thisTrcFilename,0,1);
        channelNames=output.Names;
        thisTriggerStr=triggerStrings{r};
        tmp = strcmp(channelNames, thisTriggerStr);
        indxTrigger=find(tmp);
        numfound=numel(indxTrigger);
        thisLeadField=computeLeadField(thisTrcFilename,indxTrigger,locFilename,trc2elc,threshold,fhigh,show);
        leadFields(:,r)=thisLeadField;
        
    catch
        badElecIndx=[badElecIndx r];
    end
end

%%
figure
eloc = readlocs( locFilename );
eloc_no_ref=eloc(2:end);
for i=1:size(leadFields,2)
    subplot(6,9,i);
    topoplot(leadFields(:,i),eloc_no_ref,'electrodes','off','numcontour',0);
    title(recordingNames{i})
    %pause
end

cm=jmaColors('usa');
colormap(cm);
hcb=colorbar('east');
%%
set(get(hcb,'xlabel'),'string','\muV');
cbPos=get(hcb,'Position');
set(hcb,'Position',[cbPos(1)+0.1 cbPos(2) 1.2*cbPos(3) 1.2*cbPos(4)]);

% titles={'B2-B3','B6-B7','B12-B13','A2-A3','A6-A7','A12-A13'};
% for s=1:6
%     muScalpReordered = computeLeadField(trcFilenames{s},triggerChannels(s),locFilename,trc2elc,threshold,fhigh,show);
%     subplot(2,3,s);
%     topoplot(muScalpReordered,eloc_no_ref,'electrodes','numbers');
%     title(titles{s});
% end

%%
print -dpng ../figures/leadFieldsPretty

