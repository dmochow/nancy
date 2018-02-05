function [psth,binCenters] = computePsth(mua,binWidthMs,fs)
nSamples=size(mua,1); % time dimension 1
nTrials=size(mua,2); 
binWidthSamples=binWidthMs/1000*fs;
if nSamples==binWidthSamples
    nBins=1;
else
    nBins=floor(nSamples/binWidthSamples);
end

binEdgesLeft=(0:nBins-1)*binWidthSamples+1;
binEdgesRight=(1:nBins)*binWidthSamples;
% for now ignoring any times outside of nBins*binWidthSamples (truncated
% bins)

psth=zeros(nBins,1);
for b=1:nBins
    %b/nBins
    tmp=sum(mua(binEdgesLeft(b):binEdgesRight(b),:),1);
    if nTrials>1
        psth(b)=mean(tmp);
    else
        psth(b)=tmp;
    end
end

binCenters=mean([binEdgesLeft.' binEdgesRight.'],2)/fs;