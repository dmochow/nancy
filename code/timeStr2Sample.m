function [secsElapsed,samplesElapsed] = timeStr2Sample(timeStr,refTimeStr,fs)
%   compute time elapsed between string time and reference time
if nargin<3, fs=1; end

[~, ~, ~, H, M, S] = datevec(timeStr);
secs=H*3600+M*60+S;
%samples=round(secs*fs);

[~, ~, ~, H, M, S] = datevec(refTimeStr);
refSecs=H*3600+M*60+S;
%refSamples=round(secs*fs);

secsElapsed=secs-refSecs;
samplesElapsed=round((secs-refSecs)*fs);

end

