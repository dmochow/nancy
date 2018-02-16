function [dataBlocks]  = trc2Blocks(trcfilename, blockStartTimesStr, blockLength)

%--------------------------------------------------------------------------
% trcc2Blocks reads a .TRC file and converts it into blocks defined by the
% parameters.
%
% input: 
%       trcfilename - path to .trc file
%       blockStartTimesStr - Cell array of strings that contain block start
%                            times in hh:mm:ss.ms format
%       blockLength - Length of a block in seconds
%                   
% output: 
%       dataBlocks: Blocks from the trc file based on the parameters given
%       in the form of a 3d matrix: dimensions of electrode (178), time (samples), and blocks 
%
%--------------------------------------------------------------------------
% author: Prakhyat Singh
%--------------------------------------------------------------------------

% Parameters for read_micromed function
readstim = 0;
inmv = 1;
% Read out the data from the trc file
[data, output]  = icem_read_micromed_trc(trcfilename, readstim, inmv);


numBlocks = length(blockStartTimesStr);
sampleRate = output.SR;
startTime = output.Start;
numElectrodes = output.NbrCh;
blockLengthSampleSize = blockLength * sampleRate; % Number of samples
data = data{1}; % Going from struct to a 2-D array embedded within it.

% Convert startTime to seconds
% Put hours, minutes, seconds in an array
timeArray = str2num(strrep(startTime,':',' ')); 
% Convert hours, minutes, seconds to seconds
startTimeInSeconds = timeArray(1)*3600 + timeArray(2)*60 + timeArray(3);

% Array to store the start locations of blocks present in the trc data
blockStartLocs = zeros(1, numBlocks);

for i = 1:numBlocks
    timeArray = str2num(strrep(blockStartTimesStr{i},':',' '));
    % Convert to seconds and subtract start time (start of the trace)
    blockStartLocs(i) = ((timeArray(1)*3600 + timeArray(2)*60 + timeArray(3))-startTimeInSeconds)*sampleRate;    
end

% Reshape trc data to Blocks in the form of a 3d matrix: dimensions of electrode (178), time (samples), and blocks
dataBlocks = zeros(numElectrodes, blockLengthSampleSize, numBlocks);
for i = 1:numBlocks
    dataBlocks(:,:,i) = data(:,(blockStartLocs(i):blockStartLocs(i) + blockLengthSampleSize-1));
end
