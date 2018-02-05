
function [spikes,waveform]=removeBadSpikes(spikesIn,waveformIn,thresh)
% remove spikes that don't look like spikes
if nargin<3, thresh=0.5; end
validIndx=3:58; %HARD-CODED (for now... and ever)
midIndx=25:40;  %HARD-CODED
nlenergy=nle(waveformIn);

if isempty(spikesIn)
    spikes=[];
    waveform=[];
else
    props=sum( nlenergy(:,midIndx).^2 , 2 ) ./ sum( nlenergy(:,validIndx).^2 , 2 );
    goodSpikeIndx=find(props>thresh);
    spikes=spikesIn(goodSpikeIndx,:)
    waveform=waveformIn(goodSpikeIndx,:);
end

