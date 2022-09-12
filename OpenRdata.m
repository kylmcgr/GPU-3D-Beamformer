%% Pieces of code used open data

% Open RData matrix   
% dir_save = 'D:\SURF 2022 Data\Power_Doppler_Direct_4angles_22V\' ;
% load([dir_save 'P.mat']);
% load([dir_save 'Receive.mat']);
% load([dir_save 'Trans.mat']);
% load([dir_save 'TX.mat']);
% clear RData
% for Dop_count = 1
%     for frame_num = 1
%         fid = fopen([dir_save ,sprintf('RData_Dop%.3d_frame%.3d.bin',Dop_count,frame_num)], 'r');
%         RData_frame = fread(fid, 'int16');
%         fclose(fid);
% 
%         RData_frame = reshape(RData_frame,[],256);
%         RData(Dop_count,frame_num,:,:) = int16(RData_frame);
%     end
% end


% Dop_count = 1;
% frame_num = 1;
% dir_save = '\' ;

fid = fopen('D:\SURF 2022 Data\Direct_data_6_PlaneWaves_asymetrical_pyramid\RData_Dop001_frame001.bin', 'r');
RData = fread(fid, 'int16');
fclose(fid);

RData = reshape(RData,[],256);
RData = int16(RData);



