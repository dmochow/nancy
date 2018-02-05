%clear all; close all; clc

load('../precomputed/ssvep_jj','ssvep');
%
ssvep6=ssvep(:,[1 3 5]);
ssvep1=ssvep(:,[2 4 6]);


figure
subplot(211);
plot([ssvep1(:,1) ssvep1(:,2) ssvep1(:,3)]); 
legend('pre','stim','post'); legend boxoff; title('1.2 Hz')
subplot(212);
plot([ssvep6(:,1) ssvep6(:,2) ssvep6(:,3)]);
legend('pre','stim','post'); legend boxoff; title('6 Hz')

% ssvep6_ds=(ssvep6(:,2)-ssvep6(:,1))./ssvep6(:,1);
% ssvep6_dp=(ssvep6(:,3)-ssvep6(:,1))./ssvep6(:,1);
% ssvep1_ds=(ssvep1(:,2)-ssvep1(:,1))./ssvep1(:,1);
% ssvep1_dp=(ssvep1(:,3)-ssvep1(:,1))./ssvep1(:,1);
% 
% figure;
% subplot(221); stem(ssvep6_ds);
% subplot(222); stem(ssvep6_dp);
% subplot(223); stem(ssvep1_ds);
% subplot(224); stem(ssvep1_dp);
