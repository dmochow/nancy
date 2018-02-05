function [brainDataWindow,channelNames,fs]=readStimulationData(trcFilename,fhigh,show)

if nargin<3, show=0; end
if nargin<2, fhigh=1;  end

% read the data into matlab format
readstim=0;
inmv=1;
[data,output]  = icem_read_micromed_trc(trcFilename, readstim, inmv);
fs=output.SR;
data=data{1};

% design filter (nothing below fhigh Hz + notch at 50 Hz)
[notchnum,notchdenom]=butter(2,[48 52]/fs*2,'stop');
if fhigh % anything but zero
    [hpnum,hpdenom]=butter(2,fhigh/fs*2,'high'); % high-pass filter
else
    hpnum=1; hpdenom=1;
end
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
tmp = strfind(channelNames, 'MKR');
indxMarker=find(~cellfun(@isempty,tmp));
tmp = strfind(channelNames, 'PULS');
indxPulse=find(~cellfun(@isempty,tmp));
tmp = strfind(channelNames, 'BEAT');
indxBeat=find(~cellfun(@isempty,tmp));
tmp = strfind(channelNames, 'ECG');
indxECG=find(~cellfun(@isempty,tmp));
tmp = strfind(channelNames, 'SpO2');
indxSpO2=find(~cellfun(@isempty,tmp));
indxKeep=setdiff(1:numel(channelNames),[indxScalp;indxMarker;indxPulse;indxBeat;indxECG;indxSpO2]);
channelNames=channelNames(indxKeep);
%%
brainData=data(indxKeep,:);
figure;
plot(brainData.');
title('Click on two points to indicate start and end of analysis window.');
[X,Y]=ginput(2);

if X(1)<0 & X(2)<0
    brainDataWindow=brainData; % if we don't want to window
else
    brainDataWindow=brainData(:,round(X(1)):round(X(2)));
end

if show
    figure;
    plot(brainDataWindow.');
end
%%

