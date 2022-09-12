%% Pieces of code used open IQData


% Path = 'C:\Users\User\OneDrive - California Institute of Technology\Caltech\Junior\SURF\Beamformer_scripts\Raw_data_simulation_oneplanewave_onetransmit_onereceive\';

%% SVD loop
% cd(Path)
load('P.mat')
load('Trans.mat')

% input parameters
P.Lat_dim = P.PData.Size(1) ;
P.Ax_dim = P.PData.Size(3); 
ensemble_length = P.numIQperDopplerIm;
fid = fopen('IQ3D_block_001_frame001.bin');
tmp = fread(fid, 'double'); %fread
fclose(fid);
IQ_temp = reshape(tmp,P.Lat_dim, P.Lat_dim, 2*P.Ax_dim ,ensemble_length);
IQData = IQ_temp(:,:,1:P.Ax_dim,:)+1i*IQ_temp(:,:,(P.Ax_dim)+1:P.Ax_dim*2,:);

figure
for x=1:82
    imagesc(squeeze(abs(IQData(:,:,x))))
    pause(0.1)
end

% 
% ii = 1 ;
% fid = fopen('IQ3D_block_001_frame001.bin');
% tmp = fread(fid, 'double'); %fread
% fclose(fid);
% IQ_temp = reshape(tmp,P.Lat_dim, P.Lat_dim, 2*P.Ax_dim ,ensemble_length);
% IQData = IQ_temp(:,:,1:P.Ax_dim,:)+1i*IQ_temp(:,:,(P.Ax_dim)+1:P.Ax_dim*2,:);