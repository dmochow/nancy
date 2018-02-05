function [muScalpReordered,nStimFound] = computeLeadField(trcFilename,TRIGGER_CHANNEL,locFilename,trc2elc,THRESHOLD,fhigh,show)

if nargin<7, show=0; end;
if nargin<6, fhigh=10; end
if nargin<5, THRESHOLD=2000; end 
if nargin<4
    trc2elc=[11 12 13 14 1 2 3 4 5 6 ...
    7 8 9 10 23 24 20 21 22 15 16 17 18 19];
end
if nargin<3, locFilename='../data/PAT_1/HUE_LI_EEG-SEEG2.elc'; end
if nargin<2, error('At least two argument required'); end

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

% apply filter
prependLen=round(5*fs);
data=cat(2,zeros(size(data,1),prependLen),data); % prepend zeros
data=filter(b,a,data,[],2); % apply preproc filter
data=data(:,prependLen+1:end); % undo prepend


%% now get subsets of data channels
channelNames=output.Names;
tmp = strfind(channelNames, 's');
indxScalp=find(~cellfun(@isempty,tmp));

% create scalp data
scalpData=data(indxScalp,:);
scalpData=scalpData-repmat(mean(scalpData,1),size(scalpData,1),1);% common average 

%%
trigger=data(TRIGGER_CHANNEL,:);
highTimes=find(trigger>THRESHOLD);
% while isempty(highTimes)
%     THRESHOLD=THRESHOLD*0.8;
%     highTimes=find(trigger>THRESHOLD);
% end
diffHighTimes=[0 diff(highTimes)];
stimTimes(1)=highTimes(1);
stimTimes=cat(2,stimTimes,highTimes(diffHighTimes>1));
nStimFound=numel(stimTimes);
fprintf('Found %d stimulation onset times \n',numel(stimTimes));
stimTimes=stimTimes+0; % add latency here
if abs ( mean(diff(stimTimes)) - fs   ) >  0.01*fs % throw a warning if off by more than 1 percent
    warning('Inferred stimulation times seem to be wrong');
end

if show
    ch2show=8;
    figure;
    subplot(211)
    plot(trigger);hold on
    stem(stimTimes,trigger(stimTimes))
    subplot(212)
    plot(scalpData(ch2show,:)); hold on
    stem(stimTimes,scalpData(ch2show,stimTimes))
end

%% compute mean scalp
%muScalp=mean(scalpData(:,stimTimes),2);
muScalp=trimmean(scalpData(:,stimTimes),5,2);
muScalpReordered=muScalp(trc2elc);

% show on scalp map
if show
    figure;
    eloc = readlocs( locFilename );
    eloc_no_ref=eloc(2:end);
    topoplot(muScalpReordered, eloc_no_ref,'electrodes','numbers');
end

