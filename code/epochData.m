function epochedData = epochData(data,timeRecordStart,timeBlockStart,blockDurationSeconds,fs)

blockDurationSamples=round(blockDurationSeconds*fs);

[~, ~, ~, H0, MN0, S0] = datevec(timeRecordStart);
timeStartSecs=H0*3600+MN0*60+S0;

[~, ~, ~, H1, MN1, S1] = datevec(timeBlockStart);
timeBlockStartSecs=H1*3600+MN1*60+S1;

% offset everything to onset
timeBlockStartSecs=timeBlockStartSecs-timeStartSecs;
timeBlockStartSamples=round(timeBlockStartSecs*fs);

epochedData=data(:,timeBlockStartSamples+1:timeBlockStartSamples+blockDurationSamples);
