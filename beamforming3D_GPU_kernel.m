%% GPU beamformer kernel for 3D IQ data single image
% Beamforms a single image given RData,Receive,TX,P,Trans objects
% Author - Kyle McGraw
% Written on August 2022
% 

function bf = beamforming3D_GPU_kernel(RData,Receive,TX,Trans,P);
speedOfSound=1540; % m/s
invSpeedOfSound=1/speedOfSound; % s/m
freq=Trans.frequency*1e6; % wavelengths/s
lambda=speedOfSound/freq; % m
Senscutoff=0.6; % cutoff for ElementSens value, 0.6 is default, should correspond to 9.8 degrees off axis

xPiezo = Trans.ElementPos(1:32,1)*1e-3; % x position of 1024 elements -> 32
yPiezo = Trans.ElementPos(1:32:end,2)*1e-3; % y position of 1024 elements -> 32
% zPiezo = Trans.ElementPos(:,3)*1e-3; % z position of 1024 elements, should be all 0

X=((0:P.PData.Size(1)-1)+P.PData.Origin(1))*lambda; % Size of matrix that we are imaging
Y=(-(0:P.PData.Size(2)-1)+P.PData.Origin(2))*lambda; % Size in wavelength, starting at origin, changed to m
Z=P.startDepth_mm*1e-3:lambda:P.endDepth_mm*1e-3;

% Might need to check on this, but should work for most cases of direct and
% light transmissions
% load('event_seq_mat.mat');
if exist('event_seq_mat','var')
    % If you have the event seq mat file this should be good
    TX_indices=[event_seq_mat.tx];
else
    % Lines up the tx indices to the receive indices, can use the event_seq_mat
    % instead if you have it
    if isfield(P,'nb_light_TX_RX') && P.nb_light_TX_RX > 1
        TX_indices=repmat(1:length(TX), 1, length(Receive)/(length(TX)/4*10)); % If we have a single transmit receive
        TX_indices=repelem(TX_indices, repmat([2 3 3 2],1,length(TX_indices)/4)); % If we have 1,2 1,2,3 2,3,4 3,4 recieve for 1 2 3 4 transmit
    else
        TX_indices=repmat(1:length(TX), 1, length(Receive)/length(TX)); % If we have a single transmit receive
    end
end

% Send all the data to GPU memory
bf=gpuArray(complex(single(zeros(length(X),length(Y),length(Z)))));
grid_X=gpuArray(single(X));
grid_Y=gpuArray(single(Y));
grid_Z=gpuArray(single(Z));
piezo_X=gpuArray(single(xPiezo));
piezo_Y=gpuArray(single(yPiezo));
idxs_TX=gpuArray(int32(TX_indices));
TX_apertures=gpuArray(int32([TX.aperture]));
RCV_apertures=gpuArray(int32([Receive.aperture]));
angleDelays=gpuArray(single([TX.Delay]));
startSamples=gpuArray(int32([Receive.startSample]));
endSamples=gpuArray(int32([Receive.endSample]));
data=gpuArray(int16(RData));

g = gpuDevice;
kernel_name = 'beamforming3D_GPU_cuda';
% bf_flag = 'beamforming3D_cuda_v1';
system('nvcc beamforming3D_GPU_cuda.cu -ptx -o beamforming3D_GPU_cuda.ptx'); % compileCUDA(kernel_name)
k = parallel.gpu.CUDAKernel([kernel_name '.ptx'],[kernel_name '.cu']); % ,bf_flag);

setConstantMemory(k,...
    'invSpeedOfSound',single(invSpeedOfSound),...
    'xPixLen',int32(length(X)),...
    'yPixLen',int32(length(Y)),...
    'zPixLen',int32(length(Z)),...
    'xPiezoLen',int32(length(xPiezo)),...
    'yPiezoLen',int32(length(yPiezo)),... %     'zPiezo',int32(zPiezo),...
    'receiveLen',int32(length(TX)),...
    'freq',single(freq),...
    'pi',single(pi),...
    'sampleDim',int32(size(RData,1)),...
    'sensCutoff',single(Senscutoff),...
    'elementSens',single(Trans.ElementSens)...
    );
optimizeThreadBlocSize = 0; 

if optimizeThreadBlocSize
    disp('Searching the fastest block size...')
    bestBlock = optimizeBlockSize_iq_kernel_beamformer_allFrames_100(k, g, IQout, maskout, data, buff{1}{1}, Px, Py, Pz, BF ); 
    k.ThreadBlockSize = bestBlock;
else
%     warning('trhread block size .?')
    k.ThreadBlockSize = [16 2 8];
end

dataSize = [length(X) length(Y) length(Z)];

% set grid size
k.GridSize = ceil([dataSize(1)/k.ThreadBlockSize(1),dataSize(2)/k.ThreadBlockSize(2),dataSize(3)/k.ThreadBlockSize(3)]);

k.SharedMemorySize = 0;

bf=feval(k,bf,data,grid_X,grid_Y,grid_Z,piezo_X,piezo_Y,idxs_TX,TX_apertures,RCV_apertures,angleDelays,startSamples,endSamples);
bf=gather(bf);
bf=permute(bf,[2 1 3]); % To make same dimensions as Verasonics beamformer

wait(g)