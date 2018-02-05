function [ssvep,allP,channelNames]=readFaceData(trcFilename,fhigh)

if nargin<2, fhigh=1;  end

% read the data into matlab format
readstim=0;
inmv=1;
[data,output]  = icem_read_micromed_trc(trcFilename, readstim, inmv);
fs=output.SR;
data=data{1};

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

% get subset of non-trigger/physio channels
brainData=data(indxKeep,:);
trigger=data(indxMarker(end),:);

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
brainData=cat(2,zeros(size(brainData,1),prependLen),brainData); % prepend zeros
brainData=filter(b,a,brainData,[],2); % apply preproc filter
brainData=brainData(:,prependLen+1:end); % undo prepend


figure(1);
plot(trigger);
title('Click to the left of a block start.  Click on negative x to stop.');
c=0; % counter of blocks
while 1
    c=c+1;
    [X,Y]=ginput(1);
    if X<0, break; end
    thisSelection=round(X);
    blockStart=find( trigger(thisSelection:end) < -100 ,1); % find start of stimulus
    windowStart=thisSelection+blockStart+2*fs;
    windowEnd=windowStart+75*(1/1.2)*fs-1;
    window=windowStart:windowEnd;
    brainDataWindow=brainData(:,window);
    hold on; stem(windowStart,trigger(windowStart),'ro');
    stem(windowEnd,trigger(windowEnd),'ro');
    %blockOnsets=cat(1,blockOnsets,thisSelection);
    
    % fft
    nSamples=size(brainDataWindow,2);
    nfft=2^nextpow2(nSamples);
    freqs=(0:nfft-1)/nfft*fs;
    [~,f1]=min(abs(freqs-1.2));
    [~,f2]=min(abs(freqs-6));
    
    fftBrainDataWindow=fft(brainDataWindow,nfft,2);
    P=abs(fftBrainDataWindow)/nfft;
    
    allP(:,:,c)=P;
    %ssveps(:,:,c)=P(:,[f1 f2]);
    
end


N=size(allP,3);
meanP=squeeze(mean(allP,1));
for n=1:N
    subplot(2,2,n);
    plot(freqs,meanP(:,n)); xlim([0 30]); hold on;
    stem(freqs(f1),meanP(f1,n));
    stem(freqs(f2),meanP(f2,n));
end

ssvep=allP(:,[f1 f2],:);


% brainDataWindow=[];
% for r=1:nRegions
%     brainDataWindow=cat(2,brainDataWindow,brainData(:,X((r-1)*2+1):X(2*r) ));
% end
%else
% trigger logic here
%end






