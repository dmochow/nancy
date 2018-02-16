% This script calls trc2Block (which in turn calls read_micromed_trc) to
% ... divide the trc file data into blocks as per the passed parameters.

clear all; close all; clc;

addpath('../../jd-lab-common')

trcfilename = '../../EEG_2014.TRC';
blockLength = 70;
blockStartTimesStr={'10:23:16';'10:24:42';'10:26:06';... %pre tDCS
'10:30:24'; '10:32:13'; '10:33:43'; ... % during tDCS
'10:36:16'; '10:37:57'; '10:39:18'}; % post tDCS

Result  = trc2Blocks(trcfilename, blockStartTimesStr, blockLength);