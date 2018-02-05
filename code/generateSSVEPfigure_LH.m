clear all; close all; clc;
addpath(genpath('../../COMMON'));

viewAngle=[-6 76];
zoomFactor=1.5;
faceAlpha=0.25;
nElectrodes=184;
fs=2048;
ch2show=76; % target

matFilenames={'../data/HUE_LI/Extracted tDCS/FaceLoc1ssvep.mat', ...
    '../data/HUE_LI/Extracted tDCS/FaceLoc2ssvep.mat', ...
    '../data/HUE_LI/Extracted tDCS/tDCS-FaceLoc3ssvep.mat', ...
    '../data/HUE_LI/Extracted tDCS/FaceLoc5ssvep.mat', ...
    '../data/HUE_LI/Extracted tDCS/FaceLoc6ssvep.mat'};

nFilenames=numel(matFilenames);
ssveps={};

load('../data/HUE_LI/Extracted/montages_9electrodes_1ma/montage_A1213_9_field_and_montage.mat','Il1');
locFilename='../data/HUE_LI/HUE_LI_EEG-SEEG2.elc';

surfFilename_right='../data/HUE_LI/freesurfer/surf/rh.pial';
surfFilename_left='../data/HUE_LI/freesurfer/surf/lh.pial';

[verts_right, faces_right] = read_surf(surfFilename_right);
faces_right=faces_right+1;
nVerts_right=size(verts_right,1);
vals_right=zeros(nVerts_right,1);
%
[verts_left, faces_left] = read_surf(surfFilename_left);
faces_left=faces_left+1;
nVerts_left=size(verts_left,1);
vals_left=zeros(nVerts_left,1);


electrodeCoordinateFilename='../data/HUE_LI/Extracted/electrodeCoordinates.mat';
load(electrodeCoordinateFilename,'locs','labels'); % locs and labels now in memory
nLocs=size(locs,1);
channelLocs=NaN*ones(nElectrodes,3); %

for f=1:nFilenames
    load(matFilenames{f},'ssvep','channelNames','allP');
    
    tmpL=circshift(allP,[0 1 0]);
    tmpR=circshift(allP,[0 -1 0]);
    allPnorm=allP-(tmpL+tmpR)/2;
    nfft=size(allP,2);
    freqs=(0:nfft-1)/nfft*fs;
    [~,indx6]=min(abs(freqs-6));
    [~,indx12]=min(abs(freqs-12));
    [~,indx18]=min(abs(freqs-18));
    
    allSpectra{f}=allPnorm;
    
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
end


%%
% fill in location of target electrode (!)
aRange=65:75;
mdiff=mean(diff(channelLocs(aRange,:)));
channelLocs(76,:)=channelLocs(75,:)+mdiff;

%%
ssvepPre=cat(3,ssveps{1},ssveps{2});
ssvepStim=ssveps{3};
ssvepPost=cat(3,ssveps{4},ssveps{5});

mussvepPre_jd=mean(ssvepPre,3);
mussvepStim_jd=mean(ssvepStim,3);
mussvepPost_jd=mean(ssvepPost,3);

%%
del1=(mussvepStim_jd(:,1)-mussvepPre_jd(:,1))./mussvepPre_jd(:,1);
del2=(mussvepPost_jd(:,1)-mussvepPre_jd(:,1))./mussvepPre_jd(:,1);
del3=(mussvepStim_jd(:,2)-mussvepPre_jd(:,2))./mussvepPre_jd(:,2);
del4=(mussvepPost_jd(:,2)-mussvepPre_jd(:,2))./mussvepPre_jd(:,2);



%%
% nChannels=size(mussvepPre_jd,1);
% titleStr={'1.2 Hz','6 Hz'};
% figure
% for f=1:2
%     data=[mussvepPre_jd(:,f) mussvepStim_jd(:,f) mussvepPost_jd(:,f)];
%     sems=std(data,[],1)/sqrt(size(data,1));
%         
%     hs(3)=subplot(2,3,3); hold on
%     if f==1
%         herr=errorbar(mean(data),sems,'b','LineStyle','none');
%     else
%         herr=errorbar([1.1 2.1 3.1],mean(data),sems,'r','LineStyle','none');
%     end
% %     for j=1:3
% %         hpp(j)=plot(j,mean(data(:,j)),'o','LineStyle','none');
% %         set(hpp(j),'MarkerFaceColor',get(hp(j),'Color'));
% %         set(hpp(j),'MarkerEdgeColor',get(hp(j),'Color'));
% %     end
%     
%     htit(1)=title('All Electrodes','FontWeight','normal');
%     ylabel('Amplitude (\mu V)');
%     xlim([0.5 3.5]);
%     yl=ylim;
%     ylim([yl(1) 2.6]);
%     set(gca,'Xtick',1:3);
%     set(gca,'XTickLabel',{'Pre-tDCS','tDCS','Post-tDCS'});
%     %hlg2=legend('1.2 Hz','6 Hz'); set(hlg2,'box','off');  
%     htext=text(0.8,2,'1.2 Hz');
%     htext2=text(1,2.5,'6 Hz');
%     
% end



%%
% clean figure begins here
figure
desiredE=0.5*ones(nElectrodes,1);
desiredE(ch2show)=1;
S=30;
cm=jmaColors('coolhotcortex');
brainColor=0.5;

hs(1)=subplot(231);
scatter3(channelLocs(:,1),channelLocs(:,2),channelLocs(:,3),S,desiredE,'filled');
patch( 'Vertices',verts_right, 'Faces',faces_right,...
    'FaceVertexCData',brainColor*ones(nVerts_right,1),...  %gCortex.colors,...   %
    'FaceColor','interp', 'FaceLighting','gouraud', 'BackFaceLighting','unlit', ...
    'EdgeColor','none', 'DiffuseStrength',0.7, 'SpecularStrength',0.05,...
    'SpecularExponent',5, 'SpecularColorReflectance',0.5 , ...
    'LineStyle','none','FaceAlpha',faceAlpha);
patch( 'Vertices',verts_left, 'Faces',faces_left,...
    'FaceVertexCData',brainColor*ones(nVerts_left,1),...  %gCortex.colors,...   %
    'FaceColor','interp', 'FaceLighting','gouraud', 'BackFaceLighting','unlit', ...
    'EdgeColor','none', 'DiffuseStrength',0.7, 'SpecularStrength',0.05,...
    'SpecularExponent',5, 'SpecularColorReflectance',0.5 , ...
    'LineStyle','none','FaceAlpha',faceAlpha);
caxis([0 1]);
colormap(cm);
axis equal
axis off
view([91 7]);
htit(1)=title('Targeted Electrode','FontWeight','normal');
titPos=get(htit(1),'Position');
set(htit(1),'Position',[titPos(1)+45 titPos(2) titPos(3)]);



tmpPre=cat(3,allSpectra{1},allSpectra{2});
tmpStim=allSpectra{3};
tmpPost=cat(3,allSpectra{4},allSpectra{5});
spectrumPre=mean(tmpPre(ch2show,:,:),3);
spectrumStim=mean(tmpStim(ch2show,:,:),3);
spectrumPost=mean(tmpPost(ch2show,:,:),3);

spectrumsPre=mean(tmpPre(setdiff(1:nElectrodes,ch2show),:,:),3);
spectrumsStim=mean(tmpStim(setdiff(1:nElectrodes,ch2show),:,:),3);
spectrumsPost=mean(tmpPost(setdiff(1:nElectrodes,ch2show),:,:),3);
muSpectrumsPre=mean(spectrumsPre,1);
muSpectrumsStim=mean(spectrumsStim,1);
muSpectrumsPost=mean(spectrumsPost,1);
semsSpectrumsPre=std(spectrumsPre,[],1)/sqrt(size(spectrumsPre,1));
semsSpectrumsStim=std(spectrumsStim,[],1)/sqrt(size(spectrumsPre,1));
semsSpectrumsPost=std(spectrumsPost,[],1)/sqrt(size(spectrumsPre,1));


[~,f1]=min(abs(freqs-1.2));
[~,f2]=min(abs(freqs-6));

preVals=spectrumPre([f1 f2]);
stimVals=spectrumStim([f1 f2]);
postVals=spectrumPost([f1 f2]);

hs(2)=subplot(232);
bar([preVals' stimVals' postVals'].');
ylabel('Power Spectrum (\muV / Hz)')
set(gca,'XTickLabel',{'Pre','tDCS','Post'});
hlg=legend('1.2 Hz','6 Hz'); set(hlg,'box','off');
title('Targeted Electrode','FontWeight','normal');
yl=ylim;

hs(3)=subplot(233); hold on
preVals=muSpectrumsPre([f1 f2]);
stimVals=muSpectrumsStim([f1 f2]);
postVals=muSpectrumsPost([f1 f2]);

preSems=semsSpectrumsPre([f1 f2]);
stimSems=semsSpectrumsStim([f1 f2]);
postSems=semsSpectrumsPost([f1 f2]);

bar([preVals' stimVals' postVals'].');
eps=0.15;
errorbar([1-eps 2-eps 3-eps],[preVals(1) stimVals(1) postVals(1)],[preSems(1) stimSems(1) postSems(1)],'.k');
errorbar([1+eps 2+eps 3+eps],[preVals(2) stimVals(2) postVals(2)],[preSems(2) stimSems(2) postSems(2)],'.k');
ylabel('Power Spectrum (\muV / Hz)')
set(gca,'Xtick',[1 2 3]);
set(gca,'XTickLabel',{'Pre','tDCS','Post'});
hlg=legend('1.2 Hz','6 Hz'); set(hlg,'box','off');
title('Non-Targeted Electrodes','FontWeight','normal');
ylim(yl);

print -dpng ../figures/ssvepFigureHueLiClean_A1213
%%
% 
% 
% % hs(2)=subplot(232);
% % hold on
% % plot(freqs,spectrumPre,'b'); 
% % xlim([0 15]);
% % set(gca,'XTick',[1.2 6 12 18]);
% % plot(freqs,spectrumStim,'g');
% % plot(freqs,spectrumPost,'r');
% % hlg=legend('Pre','tDCS','Post');
% % ylabel('Magnitude Spectrum (\muV)');
% % xlabel('Frequency (Hz)');
% % htit(2)=title('Targeted Electrode','FontWeight','normal');
% % set(hlg,'box','off');
% % yl=ylim;
% % ylim([0 yl(2)]);
% 
% 
% %%  DCS
% eloc = readlocs( locFilename );
% eloc_no_ref=eloc(2:end);
% hs(1)=subplot(231);
% topoplot_dc(Il1*1e3,eloc_no_ref,'electrodes','off');
% cm = jmaColors('usa');
% colormap(gca,cm); 
% hcb=colorbar('south');
% set(get(hcb,'ylabel'),'String','mA');
% cbPos=get(hcb,'Position');
% set(hcb,'Position',[cbPos(1)+0.05 cbPos(2)-0.075 0.5*cbPos(3) 0.5*cbPos(4)]);
% htit(1)=title('tDCS Montage','FontWeight','normal');
% %%
% 
% % pos=get(hs(1,1),'Position');
% % set(hs(1,1),'Position',[pos(1) pos(2) 0.5 pos(4)]);
% % 
% % pos=get(hs(1,2),'Position');
% % set(hs(1,2),'Position',[0.7 pos(2) 0.2 pos(4)]);
% % 
% % pos=get(hs(2,1),'Position');
% % set(hs(2,1),'Position',[pos(1) pos(2) 0.5 pos(4)]);
% % 
% % pos=get(hs(2,2),'Position');
% % set(hs(2,2),'Position',[0.7 pos(2) 0.2 pos(4)]);
% 
% %
% 

% 
% 
% 
% %%
% % this for trying to visualize in 3d
% load ../data/HUE_LI/redgreyblue
% cmrgb=cm;  
% 
% S=40;
% hs(4)=subplot(223); hold on
% scatter3(channelLocs(:,1),channelLocs(:,2),channelLocs(:,3),S,del1*100,'filled') 
% caxis([-100 100]);
% colormap(cmrgb);
% axis equal
% htit(4)=title('1.2 Hz','FontWeight','normal')
% axis equal
% axis off
% fCortexRight(1) = patch( 'Vertices',verts_right, 'Faces',faces_right,...
%     'FaceVertexCData',vals_right,...  %gCortex.colors,...   %
%     'FaceColor','interp', 'FaceLighting','gouraud', 'BackFaceLighting','unlit', ...
%     'EdgeColor','none', 'DiffuseStrength',0.7, 'SpecularStrength',0.05,...
%     'SpecularExponent',5, 'SpecularColorReflectance',0.5 , ...
%     'LineStyle','none','FaceAlpha',faceAlpha);
% fCortexLeft(1) = patch( 'Vertices',verts_left, 'Faces',faces_left,...
%     'FaceVertexCData',vals_left,...  %gCortex.colors,...   %
%     'FaceColor','interp', 'FaceLighting','gouraud', 'BackFaceLighting','unlit', ...
%     'EdgeColor','none', 'DiffuseStrength',0.7, 'SpecularStrength',0.05,...
%     'SpecularExponent',5, 'SpecularColorReflectance',0.5 , ...
%     'LineStyle','none','FaceAlpha',faceAlpha);
% %view(viewAngle);
% view([86 26]);
% %zoom(zoomFactor);
% hcb=colorbar('south');
% set(get(hcb,'ylabel'),'String','% change');
% cbPos=get(hcb,'Position');
% set(hcb,'Position',[cbPos(1)+0.2 cbPos(2)-0.1 cbPos(3) cbPos(4)]);
% 
% hs(5)=subplot(224); hold on
% fCortexRight(2) = patch( 'Vertices',verts_right, 'Faces',faces_right,...
%     'FaceVertexCData',vals_right,...  %gCortex.colors,...   %
%     'FaceColor','interp', 'FaceLighting','gouraud', 'BackFaceLighting','unlit', ...
%     'EdgeColor','none', 'DiffuseStrength',0.7, 'SpecularStrength',0.05,...
%     'SpecularExponent',5, 'SpecularColorReflectance',0.5 , ...
%     'LineStyle','none','FaceAlpha',faceAlpha);
% fCortexLeft(2) = patch( 'Vertices',verts_left, 'Faces',faces_left,...
%     'FaceVertexCData',vals_left,...  %gCortex.colors,...   %
%     'FaceColor','interp', 'FaceLighting','gouraud', 'BackFaceLighting','unlit', ...
%     'EdgeColor','none', 'DiffuseStrength',0.7, 'SpecularStrength',0.05,...
%     'SpecularExponent',5, 'SpecularColorReflectance',0.5 , ...
%     'LineStyle','none','FaceAlpha',faceAlpha);
% scatter3(channelLocs(:,1),channelLocs(:,2),channelLocs(:,3),S,del3*100,'filled') 
% caxis([-100 100]);
% colormap(cmrgb);
% axis equal
% axis off
% htit(5)=title('6 Hz','FontWeight','normal');
% %view(viewAngle);
% view([86 26]);
% %zoom(zoomFactor);
% 
% titPos=get(htit(4),'Position');
% set(htit(4),'Position',[titPos(1) titPos(2) titPos(3)-50]);
% titPos=get(htit(5),'Position');
% set(htit(5),'Position',[titPos(1) titPos(2) titPos(3)-50]);
% 
% %sublabel(hs,-15,-20,'FontWeight','Bold','FontSize',16);
% %print -dpng ../figures/facesEffect3D
