%% GPU beamformer kernel to open and beamform 3D IQ data
% Opens all the data and beamforms images given directory name
% Author - Kyle McGraw
% Written on August 2022
% 

function IQs = beamforming3D_GPU_doppler_pipeline(dir_save);
% Open RData matrix   
load([dir_save 'P.mat']);
load([dir_save 'Receive.mat']);
load([dir_save 'Trans.mat']);
load([dir_save 'TX.mat']);

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
if isfile([dir_save 'event_seq_mat.mat'])
    % If you have the event seq mat file this should be good
    load([dir_save 'event_seq_mat.mat']);
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

% Sets up variables to divide data into that need for each image
indices_per_IQ = length(Receive)/P.IQperRcvFrame;
rcv_apertures = [Receive.aperture];
rcv_startSamples = [Receive.startSample];
rcv_endSamples = [Receive.endSample];

g = gpuDevice;
kernel_name = 'beamforming3D_cuda_v4';
% bf_flag = 'beamforming3D_cuda_v1';
system('nvcc beamforming3D_cuda_v4.cu -ptx -o beamforming3D_cuda_v4.ptx'); % compileCUDA(kernel_name)

% fprintf('\n\t%d dopplers, %d IQ per dop, %d frame per dop, %d IQ per frame\n', P.numDopImages, P.numIQperDopplerIm, P.numIQperDopplerIm/P.IQperRcvFrame, P.IQperRcvFrame);
IQs = zeros(P.numDopImages,P.numIQperDopplerIm,length(X),length(Y),length(Z));
for dop_num = 1%:P.numDopImages
    for frame_num = 1%:P.numIQperDopplerIm/P.IQperRcvFrame
        clear RData
        fid = fopen([dir_save, sprintf('RData_Dop%.3d_frame%.3d.bin',dop_num,frame_num)], 'r');
        RData_frame = fread(fid, 'int16');
        fclose(fid);

        RData_frame = reshape(RData_frame,[],256);
        RData = int16(RData_frame);
        for rcv = 1:P.IQperRcvFrame
%             fprintf('\n\tdop: %d\n\tframe: %d\n\trcv: %d\n', dop_num, frame_num, rcv);

            % Send all the data for the image to GPU memory
            bf=gpuArray(complex(single(zeros(length(X),length(Y),length(Z)))));
            grid_X=gpuArray(single(X));
            grid_Y=gpuArray(single(Y));
            grid_Z=gpuArray(single(Z));
            piezo_X=gpuArray(single(xPiezo));
            piezo_Y=gpuArray(single(yPiezo));
            idxs_TX=gpuArray(int32(TX_indices(indices_per_IQ*(rcv-1)+1:indices_per_IQ*rcv)));
            TX_apertures=gpuArray(int32([TX.aperture]));
            RCV_apertures=gpuArray(int32(rcv_apertures(indices_per_IQ*(rcv-1)+1:indices_per_IQ*rcv)));
            angleDelays=gpuArray(single([TX.Delay]));
            startSamples=gpuArray(int32(rcv_startSamples(indices_per_IQ*(rcv-1)+1:indices_per_IQ*rcv)));
            endSamples=gpuArray(int32(rcv_endSamples(indices_per_IQ*(rcv-1)+1:indices_per_IQ*rcv)));
            data=gpuArray(int16(RData));

            k = parallel.gpu.CUDAKernel([kernel_name '.ptx'],[kernel_name '.cu']); % ,bf_flag);

            setConstantMemory(k,...
                'invSpeedOfSound',single(invSpeedOfSound),...
                'xPixLen',int32(length(X)),...
                'yPixLen',int32(length(Y)),...
                'zPixLen',int32(length(Z)),...
                'xPiezoLen',int32(length(xPiezo)),...
                'yPiezoLen',int32(length(yPiezo)),... %     'zPiezo',int32(zPiezo),...
                'receiveLen',int32(indices_per_IQ),...
                'freq',single(freq),...
                'pi',single(pi),...
                'sampleDim',int32(size(RData,3)),...
                'sensCutoff',single(Senscutoff),...
                'elementSens',single(Trans.ElementSens)...
                );
            
            k.ThreadBlockSize = [16 2 8];
 
            dataSize = [length(X) length(Y) length(Z)];

            % set grid size
            k.GridSize = ceil([dataSize(1)/k.ThreadBlockSize(1),dataSize(2)/k.ThreadBlockSize(2),dataSize(3)/k.ThreadBlockSize(3)]);
            k.SharedMemorySize = 0;

            bf=feval(k,bf,data,grid_X,grid_Y,grid_Z,piezo_X,piezo_Y,idxs_TX,TX_apertures,RCV_apertures,angleDelays,startSamples,endSamples);
            bf=gather(bf);
            bf=permute(bf,[2 1 3]); % To make same dimensions as Verasonics beamformer
            IQs(dop_num,P.IQperRcvFrame*(frame_num-1)+rcv,:,:,:) = bf;

            wait(g)
            reset(g);
        
        end
    end
end

