[~, ~, ~, H0, MN0, S0] = datevec(recordStartTimeStr);
timeStartSecs=H0*3600+MN0*60+S0;

[~, ~, ~, H1, MN1, S1] = datevec(tdcsStartTimeStr);
timeTdcsStartSecs=H1*3600+MN1*60+S1;
timeTdcsStartSamples=round((timeTdcsStartSecs-timeStartSecs)*fs);

[~, ~, ~, H2, MN2, S2] = datevec(tdcsStopTimeStr);
timeTdcsStopSecs=H2*3600+MN2*60+S2;
timeTdcsStopSamples=round((timeTdcsStopSecs-timeStartSecs)*fs);


chIndx=4;
figure; 
hold on;
plot(data(chIndx,:));  
plot(timeTdcsStartSamples,data(chIndx,timeTdcsStartSamples),'*r')
plot(timeTdcsStopSamples,data(chIndx,timeTdcsStopSamples),'*r')

%%
binWidthMs=100;
[psth,binCenters] = computePsth(mua(:,:,1).',binWidthMs,fs);
figure; hold on
plot(binCenters,psth);
stem(timesFamSecs,ones( numel(timesFamSecs) , 1),'r');