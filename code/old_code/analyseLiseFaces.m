clear all; close all; clc;
%


addpath(genpath('/Applications/freesurfer/matlab'));
surfFilename='/Applications/freesurfer/subjects/HUE_LI/surf/rh.pial';
%%
matFilenames={'../data/Extracted tDCS/FaceLoc1ssvep.mat', ...
    '../data/Extracted tDCS/FaceLoc2ssvep.mat', ...
    '../data/Extracted tDCS/tDCS-FaceLoc3ssvep.mat', ...
    '../data/Extracted tDCS/FaceLoc5ssvep.mat', ...
    '../data/Extracted tDCS/FaceLoc6ssvep.mat'};

nFilenames=numel(matFilenames);
ssveps={};

electrodeCoordinateFilename='../data/Extracted/electrodeCoordinates.mat';
load(electrodeCoordinateFilename,'locs','labels'); % locs and labels now in memory
nLocs=size(locs,1);

nElectrodes=184;
channelLocs=NaN*ones(nElectrodes,3); %

for f=1:nFilenames
    load(matFilenames{f},'ssvep','channelNames');
    
    % discard electrodes on the left side of the head
    leftChannels=[];
    for ch=1:numel(channelNames)
        if strfind(channelNames{ch},'''');
            leftChannels=cat(1,leftChannels,ch);
        end
        
        % try to find coordinates here
        channelNameStr=strrep(channelNames{ch},' ','');
        for l=1:nLocs
            if strfind(labels{l},channelNameStr) % we got a match, bitch
               channelLocs(ch,:)=locs(l,:); 
            end
        end
        
        
    end
    rightChannels=setdiff(1:numel(channelNames),leftChannels);
    rightChannels=setdiff(1:numel(channelNames),[]); % keep all
        
    ssveps{f}=ssvep(rightChannels,:,:); % keep only the channels on the right side
    %ssveps{f}=ssvep(leftChannels,:,:); % keep only the channels on the right side
    
    % try to find the coordinates here
    %strcmp(channelNames{f},
    
    
end



%%
ssvepPre=cat(3,ssveps{1},ssveps{2});
ssvepStim=ssveps{3};
ssvepPost=cat(3,ssveps{4},ssveps{5});

mussvepPre_jd=mean(ssvepPre,3);
mussvepStim_jd=mean(ssvepStim,3);
mussvepPost_jd=mean(ssvepPost,3);

%%
% del1=(mussvepPre_jd(:,1)-mussvepStim_jd(:,1))./mussvepPre_jd(:,1);
% del2=(mussvepPre_jd(:,1)-mussvepPost_jd(:,1))./mussvepPre_jd(:,1);
% del3=(mussvepPre_jd(:,2)-mussvepStim_jd(:,2))./mussvepPre_jd(:,2);
% del4=(mussvepPre_jd(:,2)-mussvepPost_jd(:,2))./mussvepPre_jd(:,2);
% mean([del1 del2 del3 del4])

del1=(mussvepStim_jd(:,1)-mussvepPre_jd(:,1))./mussvepPre_jd(:,1);
del2=(mussvepPost_jd(:,1)-mussvepPre_jd(:,1))./mussvepPre_jd(:,1);
del3=(mussvepStim_jd(:,2)-mussvepPre_jd(:,2))./mussvepPre_jd(:,2);
del4=(mussvepPost_jd(:,2)-mussvepPre_jd(:,2))./mussvepPre_jd(:,2);
mean([del1 del2 del3 del4])
%%
nChannels=size(mussvepPre_jd,1);
titleStr={'1.2 Hz','6 Hz'};
figure
for f=1:2
    data=[mussvepPre_jd(:,f) mussvepStim_jd(:,f) mussvepPost_jd(:,f)];
    
    hs(f,1)=subplot(2,2,2*f-1); hold on
    hp=plot(data,'o','LineStyle','none');
    for i=1:3, set(hp(i),'MarkerFaceColor',get(hp(i),'Color')); set(hp(i),'MarkerSize',3); end
    title(titleStr{f});
    xlabel('Electrode index');
    ylabel('Amplitude (\mu V)');
    xlim([1 nChannels]);
    if f==1
        legend('Pre-TES','TES','Post-TES'); legend boxoff;
    end
    
    hs(f,2)=subplot(2,2,2*f); hold on
    sems=std(data,[],1)/sqrt(size(data,1));
    herr=errorbar(mean(data),sems,'k','LineStyle','none');
    for j=1:3
        hpp(j)=plot(j,mean(data(:,j)),'o','LineStyle','none');
        set(hpp(j),'MarkerFaceColor',get(hp(j),'Color'));
        set(hpp(j),'MarkerEdgeColor',get(hp(j),'Color'));
    end
    
    ht=0.01;
    base=2.55;
    hsig1=plot([1 2],[base base],'-k');
    hsig2=plot([1 1],[base-ht base+ht],'-k');
    hsig3=plot([2 2],[base-ht base+ht],'-k');
    text(1.5,base,'**');
    base=2.5;
    hsig1=plot([1 3],[base base],'-k');
    hsig2=plot([1 1],[base-ht base+ht],'-k');
    hsig3=plot([3 3],[base-ht base+ht],'-k');
    if f==1
        text(2,base,'**');
    else
        text(2,base,'*');
    end
    
    title(titleStr{f});
    ylabel('Amplitude (\mu V)');
    xlim([0.5 3.5]);
    yl=ylim;
    ylim([yl(1) 2.6]);
    set(gca,'Xtick',1:3);
    set(gca,'XTickLabel',{'Pre-TES','TES','Post-TES'});
    
end

pos=get(hs(1,1),'Position');
set(hs(1,1),'Position',[pos(1) pos(2) 0.5 pos(4)]);

pos=get(hs(1,2),'Position');
set(hs(1,2),'Position',[0.7 pos(2) 0.2 pos(4)]);

pos=get(hs(2,1),'Position');
set(hs(2,1),'Position',[pos(1) pos(2) 0.5 pos(4)]);

pos=get(hs(2,2),'Position');
set(hs(2,2),'Position',[0.7 pos(2) 0.2 pos(4)]);

sublabel([hs(1,1) hs(1,2) hs(2,1) hs(2,2)],-20,-20,'FontWeight','Bold','FontSize',16);
print -depsc ../figures/faces_LH

%%
[h,p1]=ttest(mussvepPre_jd(:,1)-mussvepStim_jd(:,1));
[h,p2]=ttest(mussvepPre_jd(:,1)-mussvepPost_jd(:,1));
[h,p3]=ttest(mussvepPre_jd(:,2)-mussvepStim_jd(:,2));
[h,p4]=ttest(mussvepPre_jd(:,2)-mussvepPost_jd(:,2));
[p1 p2 p3 p4]
%%
% this for generating the ssvep data that is used in this script
%trcFilename='../data/Extracted tDCS/FaceLoc6.trc';
% trcFilename='../data/Extracted tDCS/tDCS-FaceLoc3.trc';
% fhigh=1;
% saveFilename=[trcFilename(1:end-4) 'ssvep.mat'];
% [ssvep,allP,channelNames]=readFaceData(trcFilename,fhigh);
%save(saveFilename,'ssvep','allP','channelNames');


addpath(genpath('/Applications/freesurfer/matlab'));
surfFilename='/Applications/freesurfer/subjects/HUE_LI/surf/rh.pial';
[verts, faces] = read_surf(surfFilename);
faces=faces+1;
nVerts=size(verts,1);
vals=zeros(nVerts,1);

figure



%%
% this for trying to visualize in 3d
load ../data/redgreyblue
figure
S=40;
subplot(121); hold on
scatter3(channelLocs(:,1),channelLocs(:,2),channelLocs(:,3),S,del1*100,'filled') 
caxis([-100 100]);
colormap(cm);
axis equal
title('1.2 Hz')

fCortex = patch( 'Vertices',verts, 'Faces',faces,...
    'FaceVertexCData',vals,...  %gCortex.colors,...   %
    'FaceColor','interp', 'FaceLighting','gouraud', 'BackFaceLighting','unlit', ...
    'EdgeColor','none', 'DiffuseStrength',0.7, 'SpecularStrength',0.05,...
    'SpecularExponent',5, 'SpecularColorReflectance',0.5 , ...
    'LineStyle','none','FaceAlpha',0.25);
axis equal;
camlight('headlight')


subplot(122);
scatter3(channelLocs(:,1),channelLocs(:,2),channelLocs(:,3),S,del3*100,'filled') 
caxis([-100 100]);
colormap(cm);
axis equal
title('6 Hz');

%% 
% address laurent's comment about looking at electrodes with a strong ssvep
[sortvals_1,sortind_1]=sort(mussvepPre_jd(:,1),'descend'); % strongest ssvep at 1.2 Hz (before TDCS)
[sortvals_2,sortind_2]=sort(mussvepPre_jd(:,2),'descend'); % strongest ssvep at 6 Hz  (before TDCS)
Q=10; % how many of the strongest channels to look at
channelNames(sortind_1(1:Q))
channelNames(sortind_2(1:Q))
mean([ del1(sortind_1(1:Q)) del2(sortind_1(1:Q)) del3(sortind_2(1:Q)) del4(sortind_2(1:Q)) ])



%%
