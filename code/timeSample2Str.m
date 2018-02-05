function timeStr = timeSample2Str(refTimeStr,timeElapsedSecs,fs)
%   compute time elapsed between string time and reference time
if nargin<3, fs=1; end

[~, ~, ~, H, M, S] = datevec(refTimeStr);
refSecs=H*3600+M*60+S;
%refSamples=round(secs*fs);

timeSecs=refSecs+timeElapsedSecs;

hours=floor(timeSecs/3600);

tmp=timeSecs-hours*3600;
minutes=floor(tmp/60);

tmp2=tmp-minutes*60;

seconds=floor(tmp2);

milliseconds=(tmp2-seconds)*1000;

hstr=num2str(hours,'%0.2d');
mstr=num2str(minutes,'%0.2d');
sstr=num2str(seconds,'%0.2d');
msstr=num2str(milliseconds,'%0.3g');
timeStr=[hstr ':' mstr ':' sstr '.' msstr];
end


