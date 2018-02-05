clear all; close all; clc;

% this for generating the ssvep data
%trcFilename='../data/Extracted tDCS/FaceLoc6.trc';
% trcFilename='../data/Extracted tDCS/tDCS-FaceLoc3.trc';
% fhigh=1;
% saveFilename=[trcFilename(1:end-4) 'ssvep.mat'];
% [ssvep,allP,channelNames]=readFaceData(trcFilename,fhigh);
%save(saveFilename,'ssvep','allP','channelNames');

%%
matFilenames={'../data/Extracted tDCS/FaceLoc1ssvep.mat', ...
    '../data/Extracted tDCS/FaceLoc2ssvep.mat', ...
    '../data/Extracted tDCS/tDCS-FaceLoc3ssvep.mat', ...
    '../data/Extracted tDCS/FaceLoc5ssvep.mat', ...
    '../data/Extracted tDCS/FaceLoc6ssvep.mat'};

nFilenames=numel(matFilenames);
ssveps={};
for f=1:nFilenames
    load(matFilenames{f},'ssvep','channelNames');
    ssveps{f}=ssvep;
end

%%
ssvepPre=cat(3,ssveps{1},ssveps{2});
ssvepStim=ssveps{3};
ssvepPost=cat(3,ssveps{4},ssveps{5});

mussvepPre_jd=mean(ssvepPre,3);
mussvepStim_jd=mean(ssvepStim,3);
mussvepPost_jd=mean(ssvepPost,3);

%%
del1=(mussvepPre_jd(:,1)-mussvepStim_jd(:,1))./mussvepPre_jd(:,1);
del2=(mussvepPre_jd(:,1)-mussvepPost_jd(:,1))./mussvepPre_jd(:,1);
del3=(mussvepPre_jd(:,2)-mussvepStim_jd(:,2))./mussvepPre_jd(:,2);
del4=(mussvepPre_jd(:,2)-mussvepPost_jd(:,2))./mussvepPre_jd(:,2);
mean([del1 del2 del3 del4])

%%
nChannels=size(mussvepPre_jd,1);
figure
subplot(211);
plot([mussvepPre_jd(:,1) mussvepStim_jd(:,1) mussvepPost_jd(:,1)],'o','LineStyle','none');
% stem(1:nChannels,mussvepPre_jd(:,1)); hold on;
% stem((1:nChannels)+0.1,mussvepStim_jd(:,1));
% stem((1:nChannels)+0.2,mussvepPost_jd(:,2));
legend('pre','stim','post'); legend boxoff; title('1.2 Hz')

subplot(212);
plot([mussvepPre_jd(:,2) mussvepStim_jd(:,2) mussvepPost_jd(:,2)],'o','LineStyle','none');
legend('pre','stim','post'); legend boxoff; title('6 Hz')

%% from Jacques
load('../precomputed/ssvep_jj','ssvep');
ssvep6=ssvep(:,[1 3 5]);
ssvep1=ssvep(:,[2 4 6]);
mussvepPre=[ssvep1(:,1) ssvep6(:,1)];
mussvepStim=[ssvep1(:,2) ssvep6(:,2)];
mussvepPost=[ssvep1(:,3) ssvep6(:,3)];

del1=(mussvepPre(:,1)-mussvepStim(:,1))./mussvepPre(:,1);
del2=(mussvepPre(:,1)-mussvepPost(:,1))./mussvepPre(:,1);
del3=(mussvepPre(:,2)-mussvepStim(:,2))./mussvepPre(:,2);
del4=(mussvepPre(:,2)-mussvepPost(:,2))./mussvepPre(:,2);
mean([del1 del2 del3 del4])

figure
subplot(211);
plot([mussvepPre(:,1) mussvepStim(:,1) mussvepPost(:,1)]);
%stem([mussvepPre_jd(:,1) mussvepStim_jd(:,1) mussvepPost_jd(:,1)]);
legend('pre','stim','post'); legend boxoff; title('1.2 Hz')
subplot(212);
plot([mussvepPre(:,2) mussvepStim(:,2) mussvepPost(:,2)]);
legend('pre','stim','post'); legend boxoff; title('6 Hz')


if 0
 
%%
% reorder according to lead field matrix to facilitate comparison with
% electric field
% for each intracerebral stimulation site, get the name of the two adjacent electrodes that were stimulated
leadFieldChannelsFilename='../data/Extracted/leadFieldChannels.mat';
load(leadFieldChannelsFilename);
nLeadFieldChannels=52; %
ressvepPre=zeros(nLeadFieldChannels,2);
ressvepStim=zeros(nLeadFieldChannels,2);
ressvepPost=zeros(nLeadFieldChannels,2);
E=zeros(nLeadFieldChannels,1);
for e=1:nLeadFieldChannels
    channel1str=leadFieldChannel1{e};
    channel2str=leadFieldChannel2{e};
    tmp1 = strcmp(channelNames, channel1str);
    tmp2 = strcmp(channelNames, channel2str);
    indxChannel1=find(tmp1);
    indxChannel2=find(tmp2);
    
    ressvepPre(e,:)=0.5*(1*mussvepPre(indxChannel1,:)+1*mussvepPre(indxChannel2,:));
    ressvepStim(e,:)=0.5*(1*mussvepStim(indxChannel1,:)+1*mussvepStim(indxChannel2,:));
    ressvepPost(e,:)=0.5*(1*mussvepPost(indxChannel1,:)+1*mussvepPost(indxChannel2,:));
    
%     ressvepPre(e,:)=0.5*(mussvepPre(indxChannel1,:)-mussvepPre(indxChannel2,:));
%     ressvepStim(e,:)=0.5*(mussvepStim(indxChannel1,:)-mussvepStim(indxChannel2,:));
%     ressvepPost(e,:)=0.5*(mussvepPost(indxChannel1,:)-mussvepPost(indxChannel2,:));
    
end


%%

%%


%%
% compare with electric field
modelDataFilename='../data/Extracted/montages_9electrodes_1ma/montage_A1213_9_field_and_montage.mat';
load(modelDataFilename);

% test at f1
delEmp1=(ressvepPre(:,1)-ressvepStim(:,1));
delEmp2=(ressvepPre(:,2)-ressvepStim(:,2));
delEmp3=(ressvepPre(:,1)-ressvepPost(:,1));
delEmp4=(ressvepPre(:,2)-ressvepPost(:,2));
delMod=El1;

[r1,p1]=corrcoef(delEmp1,delMod)
[r2,p2]=corrcoef(delEmp2,delMod)
[r3,p3]=corrcoef(delEmp3,delMod)
[r4,p4]=corrcoef(delEmp4,delMod)

%%
figure;
scatter(zscore(delEmp2),zscore(delMod));

%%
figure;
subplot(211);
plot([ressvepPre(:,2)-ressvepStim(:,2)]);
subplot(212);
plot(El1);

end
