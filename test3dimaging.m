%% Pieces of code used image results

% load Data for Beamformer3DnotGPU_script
% 
% figure, for x=1:98; imagesc(squeeze(Dop(x,:,:))'), title(x), colormap hot, caxis([0 8e11]); pause(0.1), end
% figure, for x=1:98; imagesc(squeeze(bf_dop(x,:,:))'), title(x), colormap hot; pause(0.2), end

% imagesc(squeeze(abs(Dop(40,:,:)))'), colormap hot
% figure
% imagesc(squeeze(abs(bf_dop(40,:,:)))'),colormap hot

% imagesc(squeeze(abs(bf(90,:,:)))'), 
% figure
% imagesc(squeeze(abs(bf2(90,:,:)))')

figure, 
for x=1:98
    subplot(1,2,1)
    imagesc(squeeze(abs(IQData(:,x,:)))'), 
    subplot(1,2,2)
    imagesc(squeeze(abs(bf(:,x,:)))'), 
%     subplot(1,2,1)
%     imagesc(squeeze(abs(bf(x,:,:)))'), 
%     subplot(1,2,2)
%     imagesc(squeeze(abs(IQData(x,:,:)))'), 
    sgtitle(x)
%     subplot(1,2,1)
%     imagesc(squeeze(abs(bf(x,:,:)))'), 
%     subplot(1,2,2)
%     imagesc(squeeze(abs(IQData(x,:,:)))'), 
%     subplot(1,3,3)
%     imagesc(squeeze(abs(bf2(x,:,:)))')
%     subplot(3,2,3)
%     imagesc(squeeze(abs(bf(:,x,:)))'), 
%     subplot(3,2,4)
%     imagesc(squeeze(abs(IQData(:,x,:)))')
%     if x<93
%         subplot(3,2,5)
%         imagesc(squeeze(abs(bf(:,:,x)))'), 
%         subplot(3,2,6)
%         imagesc(squeeze(abs(IQData(:,:,x)))')
%     end
    pause(0.1), 
end

% figure(1)
% imagesc((b_mpw_iq(:,:,10))*20); colormap gray;
% figure(2)
% imagesc((b_mpw_iq(:,:,75))*20); colormap gray;
% figure(3)
% imagesc((b_mpw_iq(:,:,125))*20); colormap gray;
% figure(4)
% imagesc((b_mpw_iq(:,:,200))*20); colormap gray;