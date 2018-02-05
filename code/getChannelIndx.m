function channelIndx = getChannelIndx(channelNameStr,channelNum,channelNames)
% given something like 'TB1', return the index in cell array channelNames
% that contains that channel name
channelStr=cat(2,channelNameStr,num2str(channelNum));
% now loop through channel names
channelNamesNoSpace=strrep(channelNames,' ','');
isFound=strcmp(channelNamesNoSpace,channelStr);

% logic for handling results of strcmp()
if sum(isFound)==0
    channelIndx=NaN;
elseif sum(isFound)==1
    channelIndx=find(isFound,1);
else
    warning('multiple channels matching queried string found');
end

end

