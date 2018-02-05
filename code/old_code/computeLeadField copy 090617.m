clear all; close all; clc

fhigh=10; % high-pass filter cutoff frequency

trcFilename='../data/PAT_1/EEG_1.TRC';
TRIGGER_CHANNEL=81; % needs to be adapted per file, unfortunately

%trcFilename='../data/PAT_1/EEG_3.TRC';
%TRIGGER_CHANNEL=85; % needs to be adapted per file, unfortunately

%trcFilename='../data/PAT_1/EEG_5.TRC';
%TRIGGER_CHANNEL=90; % needs to be adapted per file, unfortunately

%trcFilename='../data/PAT_1/EEG_7.TRC';
%TRIGGER_CHANNEL=66; % needs to be adapted per file, unfortunately

%trcFilename='../data/PAT_1/EEG_9.TRC';
%TRIGGER_CHANNEL=70; % needs to be adapted per file, unfortunately

%trcFilename='../data/PAT_1/EEG_11.TRC';
%TRIGGER_CHANNEL=76; % needs to be adapted per file, unfortunately

trc2elc=[11 12 13 14 1 2 3 4 5 6 ...
    7 8 9 10 23 24 20 21 22 15 16 17 18 19];
%     'sFT10'   1
%     'sP10'    2    
%     'sF8'     3
%     'sT4'     4
%     'sT6'     5
%     'sFP2'    6    
%     'sF4'     7
%     'sC4'     8
%     'sP4'     9
%     'sO2'     10
%     'sFZ'     11
%     'sCZ'     12
%     'sPZ'     13
%     'sOZ'     14
%     'sFP1'    15 
%     'sF3'     16
%     'sC3'     17
%     'sP3'     18
%     'sO1'     19
%     'sF7'     20
%     'sT3'     21
%     'sT5'     22
%     'sFT9'    23    
%     'sP9'     24


% trc2elc=[6 7 8 9 10 11 12 13 14 15 2 ...
%     3 4 5 21 22 23 24 25 18 19 20 ...
%     16 17]-1;
% sFPZ	sFZ	sCz	sPZ	sOZ	sFT10 sP10	sF8	sT8	sP8	sFP2	
% sF4	sC4	sP4	sO2	sFT9 sP9 sF7	sT7	sP7	sFP1 sF3	
% sC3	sP3	sO1


THRESHOLD=2000; % in microvolts
locFilename='../data/PAT_1/HUE_LI_EEG-SEEG2.elc';

% read the data into matlab format
readstim=0;
inmv=1;
[data,output]  = icem_read_micromed_trc(trcFilename, readstim, inmv);
fs=output.SR;
data=data{1};

% design filter (nothing below 10 Hz because we only care about the high
% frequency stimulation)
[notchnum,notchdenom]=butter(2,[48 52]/fs*2,'stop');
[hpnum,hpdenom]=butter(2,fhigh/fs*2,'high'); % high-pass filter
a = poly([roots(notchdenom);roots(hpdenom)]);
b = conv(hpnum,notchnum);
%freqz(b,a,[],fs);

% apply filter
prependLen=round(5*fs);
data=cat(2,zeros(size(data,1),prependLen),data); % prepend zeros
data=filter(b,a,data,[],2); % apply preproc filter
data=data(:,prependLen+1:end); % undo prepend


%% now get subsets of data channels
channelNames=output.Names;
tmp = strfind(channelNames, 'MKR');
indxMarkers=find(~cellfun(@isempty,tmp));

tmp = strfind(channelNames, 's');
indxScalp=find(~cellfun(@isempty,tmp));
scalpChannelNames=channelNames(indxScalp);

tmp = strfind(channelNames, 'B');
indxB=find(~cellfun(@isempty,tmp));

tmp = strfind(channelNames, 'A');
indxA=find(~cellfun(@isempty,tmp));

% create scalp data
scalpData=data(indxScalp,:);
% common average (?)
scalpData=scalpData-repmat(mean(scalpData,1),size(scalpData,1),1);

%%
trigger=data(TRIGGER_CHANNEL,:);
highTimes=find(trigger>THRESHOLD);
diffHighTimes=[0 diff(highTimes)];
stimTimes(1)=highTimes(1);
stimTimes=cat(2,stimTimes,highTimes(diffHighTimes>1));
fprintf('Found %d stimulation onset times \n',numel(stimTimes));
% add latency here
stimTimes=stimTimes+0;
if abs ( mean(diff(stimTimes)) - fs   ) >  0.01*fs % throw a warning if off by more than 1 percent
    warning('Inferred stimulation times seem to be wrong');
end
ch2show=16;
%figure; plot(scalpData(ch2show,:).');
figure; 
subplot(211)
plot(trigger);hold on
stem(stimTimes,trigger(stimTimes))
subplot(212)
plot(scalpData(ch2show,:)); hold on
stem(stimTimes,scalpData(ch2show,stimTimes))


%% compute mean scalp
%muScalp=mean(scalpData(:,stimTimes),2);
muScalp=trimmean(scalpData(:,stimTimes),5,2);
muScalpReordered=muScalp(trc2elc);

% show on scalp map
figure;
eloc = readlocs( locFilename );
eloc_no_ref=eloc(2:end);
topoplot(muScalpReordered, eloc_no_ref,'electrodes','numbers');


