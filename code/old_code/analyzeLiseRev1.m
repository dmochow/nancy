clear all; close all; clc
% include all channels, and draw the e-field in three-dimensional space
dataFromDisk=1;

% for each intracerebral stimulation site, get the name of the two adjacent electrodes that were stimulated
leadFieldChannelsFilename='../data/HUE_LI/Extracted/leadFieldChannels.mat';

% the actual data file recorded during DC stimulation
trcFilenames={'../data/HUE_LI/Extracted tDCS/tDCS-SIN-G1-2-trial 2.trc'; ...
    '../data/HUE_LI/Extracted tDCS/tDCS-SIN_G9-10 - trial 2.TRC'; ...
    '../data/HUE_LI/Extracted tDCS/tDCS-SIN_G16-17.trc'; ...
    '../data/HUE_LI/Extracted tDCS/tDCS-SIN-P1-2.TRC'; ...
    '../data/HUE_LI/Extracted tDCS/tDCS-SIN-P8-9.TRC'};
indxTargets=[21;22;24;7;9];
targetsStr={'G 2','G 10','G 17','P 2','P 9'};
strTargets={'G1-G2','G9-G10','G16-G17','P1-P2','P8-P9'};
precomputedFilenames={'../precomputed/HUE_LI/tDCS-SIN-G1-2-trial 2.mat'; ...
    '../precomputed/HUE_LI/tDCS-SIN_G9-10 - trial 2.mat'; ...
    '../precomputed/HUE_LI/tDCS-SIN_G16-17.mat'; ...
    '../precomputed/HUE_LI/tDCS-SIN-P1-2.mat'; ...
    '../precomputed/HUE_LI/tDCS-SIN-P8-9.mat'};

nFiles=numel(trcFilenames);
fs=2048;
fHigh=5;
show=0;
nLeadFieldChannels=52;
nElectrodes=184;
delx=3.5;  % spacing between adjacent contacts in mm
channelLocs=NaN*ones(nElectrodes,3); %

electrodeCoordinateFilename='../data/HUE_LI/Extracted/electrodeCoordinates.mat';
load(electrodeCoordinateFilename,'locs','labels'); % locs and labels now in memory
nLocs=size(locs,1);


%% 
% draw montage first
load ../data/HUE_LI/Extracted/montages_9electrodes_1ma/montage_P89_9_field_and_montage Il1
Il1=Il1*0.5;
locFilename='../data/HUE_LI/HUE_LI_EEG-SEEG2.elc';
eloc = readlocs( locFilename );
eloc_no_ref=eloc(2:end);
figure
hs(1)=subplot(121);
topoplot_dc(Il1*1e3,eloc_no_ref,'electrodes','off','plotrad',0.7);
cm = jmaColors('usa');
colormap(gca,cm); 
hcb=colorbar('south');
set(get(hcb,'ylabel'),'String','mA');
cbPos=get(hcb,'Position');
set(hcb,'Position',[cbPos(1) cbPos(2)-0.15 cbPos(3) cbPos(4)]);
htit(1)=title('tDCS Montage','FontWeight','normal');


%%





for t=5%nFiles
    
    trcFilename=trcFilenames{t};
    indTarget=indxTargets(t);
    precomputedFilename=precomputedFilenames{t};
    
    if ~dataFromDisk % select the range yourself
        [brainDataWindow,channelNames]=readStimulationData(trcFilename,fHigh,show);
        save(precomputedFilename,'brainDataWindow','channelNames');
        % save the data here so that you don't have to keep selecting the range
    else % grab epoched data from file
        load(precomputedFilename,'brainDataWindow','channelNames')
    end
    
    durSamples=size(brainDataWindow,2);
    nfft=durSamples;
    freqs=(0:nfft-1)/nfft*fs;
    [~,indx6Hz]=min(abs(freqs-6));
    fftBrainDataWindow=fft(brainDataWindow,nfft,2);
    %Pss=abs(fftBrainDataWindow);
    %Phi_ss=angle(fftBrainDataWindow);
    
    P2=abs(fftBrainDataWindow)/durSamples;
    P1=P2(:,1:durSamples/2+1);
    P1(:,2:end-1)=2*P1(:,2:end-1);
    Pss=P1;
    
    %% convert electric potentials to electric field strength
    %load(leadFieldChannelsFilename);
    %E=zeros(nLeadFieldChannels,1);
    for e=1:nElectrodes
        
        channelName=channelNames{e};
        channelName=strrep(channelName,' ','');
        
        if isletter(channelName(1)) && isletter(channelName(2))
            if channelName(3)==''''
                channelLetter=channelName(1:2);
                channelNumStr=channelName(4:end);
            else
                channelLetter=channelName(1:2);
                channelNumStr=channelName(3:end);
            end
        else
            if channelName(2)==''''
                channelLetter=channelName(1);
                channelNumStr=channelName(3:end);
            else
                channelLetter=channelName(1);
                channelNumStr=channelName(2:end);
            end
            
        end
        
        
        if isempty(strfind(channelName,'''')) % right side of head
            channelNum=str2num(channelNumStr);
            channelNumRef=channelNum+1;
            channelNameRefNoNum=num2str(channelNumRef);
            channelNameRef=[channelLetter ' ' channelNameRefNoNum];
            tmp = strcmp(channelNames, channelNameRef);
            indxChannelRef=find(tmp,1);
            if isempty(indxChannelRef)
                % remove space and try again
                channelNameRef=strrep(channelNameRef,' ','');
                tmp = strcmp(channelNames, channelNameRef);
                indxChannelRef=find(tmp,1);
                if isempty(indxChannelRef)
                    disp(channelNames{e})
                end
            end
        else % left side of head
            channelNum=str2num(channelNumStr);
            channelNumRef=channelNum+1;
            channelNameRefNoNum=num2str(channelNumRef);
            channelNameRef=[channelLetter ''' ' channelNameRefNoNum];
            tmp = strcmp(channelNames, channelNameRef);
            indxChannelRef=find(tmp,1);
            if isempty(indxChannelRef)
                % remove space and try again
                channelNameRef=strrep(channelNameRef,' ','');
                tmp = strcmp(channelNames, channelNameRef);
                indxChannelRef=find(tmp,1);
                if isempty(indxChannelRef)
                    disp(channelNames{e})
                end
            end
        end
        
        %disp([ channelNames{e},channelNames{indxChannelRef}])

        % compute electric field using magnitude only
        if ~isempty(indxChannelRef)
            E(e)=0.001*(Pss(e,indx6Hz)-Pss(indxChannelRef,indx6Hz))/(delx);
        else
            E(e)=NaN;
        end
        
               % try to find coordinates here
        channelNameStr=strrep(channelNames{e},' ','');
        for l=1:nLocs
            if strfind(labels{l},channelNameStr) % we got a match, bitch
               channelLocs(e,:)=locs(l,:); 
            end
        end
    
        
        
        
    end
    
    indxTarget=find(strcmp(channelNames,targetsStr{t}));
    
end


%% hack to fix missing channelLocs
[D]=mean(diff(channelLocs(1:8,:),1));
channelLocs(9,:)=channelLocs(8,:)+[D(1) D(2) D(3)];

%%
surfFilename_right='../data/HUE_LI/freesurfer/surf/rh.pial';
[verts_right, faces_right] = read_surf(surfFilename_right);
faces_right=faces_right+1;
nVerts_right=size(verts_right,1);
vals_right=0.25*ones(nVerts_right,1);
%
surfFilename_left='../data/HUE_LI/freesurfer/surf/lh.pial';
[verts_left, faces_left] = read_surf(surfFilename_left);
faces_left=faces_left+1;
nVerts_left=size(verts_left,1);
vals_left=0.25*ones(nVerts_left,1);

faceAlpha=0.25;
zoomFactor=1.5;

%%
[c,d]=max((abs(E)));
if E(d)<0
    E=-E;
end
%load ../data/redgreyblue
cm=jmaColors('coolhotcortex');

hs(2)=subplot(122);
S=40;
%subplot(121);
scatter3(channelLocs(:,1),channelLocs(:,2),channelLocs(:,3),S,E*100,'filled'); 
hold on
ht=scatter3(channelLocs(indxTarget,1),channelLocs(indxTarget,2),channelLocs(indxTarget,3),2*S,E(indxTarget)*100,'filled','d'); 
caxis([min(E) max(E)]);
colormap(cm);
axis equal
axis off


fCortexRight = patch( 'Vertices',verts_right, 'Faces',faces_right,...
    'FaceVertexCData',vals_right,...  %gCortex.colors,...   %
    'FaceColor','interp', 'FaceLighting','gouraud', 'BackFaceLighting','unlit', ...
    'EdgeColor','none', 'DiffuseStrength',0.7, 'SpecularStrength',0.05,...
    'SpecularExponent',5, 'SpecularColorReflectance',0.5 , ...
    'LineStyle','none','FaceAlpha',faceAlpha);
fCortexLeft = patch( 'Vertices',verts_left, 'Faces',faces_left,...
    'FaceVertexCData',vals_left,...  %gCortex.colors,...   %
    'FaceColor','interp', 'FaceLighting','gouraud', 'BackFaceLighting','unlit', ...
    'EdgeColor','none', 'DiffuseStrength',0.7, 'SpecularStrength',0.05,...
    'SpecularExponent',5, 'SpecularColorReflectance',0.5 , ...
    'LineStyle','none','FaceAlpha',faceAlpha);
view([86 26]);
zoom(zoomFactor);
hcb=colorbar('south');
set(get(hcb,'ylabel'),'String','V/m');
cbPos=get(hcb,'Position');
set(hcb,'Position',[cbPos(1) cbPos(2)+0.025 cbPos(3) cbPos(4)]);
htit(2)=title('Intracerebral Electric Field','FontWeight','normal');
titPos=get(htit(2),'Position');
set(htit(2),'Position',[titPos(1)+45 titPos(2) titPos(3)]);

htxt=text(48.67,13.94,-15.29,'target');

sublabel(hs,45,0,'FontWeight','bold');

% filename='../figures/lise_P8P9.png';
% print('-dpng',filename);
% crop(filename);