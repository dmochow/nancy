function [letterStr,channelNum] = extractChannelLetter(channelStr)
% given something like A1 or TB1, extract 'A' or 'TB'
% grab the channel number while we're at it
channelStr=strrep(channelStr,' ','');
letterStr=channelStr(isletter(channelStr) | channelStr=='''');
channelNum=str2num(channelStr(~ (isletter(channelStr)| channelStr=='''') ));


