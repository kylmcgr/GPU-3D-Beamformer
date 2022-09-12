%% CPU beamformer for 3D IQ data
% Beamforms image given RData,Receive,TX,P,Trans objects
% Author - Kyle McGraw
% Written on August 2022
% 

function bf=beamforming3D(RData,Receive,TX,P,Trans);
speedOfSound=1540; % m/s
invSpeedOfSound=1/speedOfSound; % s/m
freq=Trans.frequency*1e6; % wavelengths/s
lambda=speedOfSound/freq; % m
Senscutoff=0.6; % cutoff for ElementSens value, 0.6 is default, should correspond to 9.8 degrees off axis

xPiezo = Trans.ElementPos(:,1)*1e-3; % x position of 1024 elements
yPiezo = Trans.ElementPos(:,2)*1e-3; % y position of 1024 elements
zPiezo = Trans.ElementPos(:,3)*1e-3; % z position of 1024 elements, should be all 0

X=((0:P.PData.Size(1)-1)+P.PData.Origin(1))*lambda; % Size of matrix that we are imaging
Y=(-(0:P.PData.Size(2)-1)+P.PData.Origin(2))*lambda; % Size in wavelength, starting at origin, changed to m
Z=P.startDepth_mm*1e-3:lambda:P.endDepth_mm*1e-3;

bf=zeros(length(X),length(Y),length(Z),'single'); % Final image matrix

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

for x_pix=1:length(X) % Loop through all pixels in imaging matrix
%     disp(x_pix);
    dX=(X(x_pix)-xPiezo)*invSpeedOfSound; % X delay from piezo to pixel 
    for y_pix=1:length(Y)
%         disp(y_pix);
        dY=(Y(y_pix)-yPiezo)*invSpeedOfSound; % Y delay from piezo to pixel
        for z_pix=1:length(Z)
            dZ=(Z(z_pix)-zPiezo)*invSpeedOfSound; % Z delay from piezo to pixel
            Delay=sqrt(dX.^2+dY.^2+dZ.^2)'; % Get 1024 array of delays from each piezo to the pixel
            elementAngles=atan(sqrt(dX.^2+dY.^2)./dZ);
            for idx=1:length(TX) % For each transmit-receive step
                TX_idx=TX_indices(idx); % Index in TX object
                tx_apod=logical(TX(TX_idx).Apod); % Transmit apod
                rcv_apod=logical(Receive(idx).Apod); % Receive apod
                angleDelayFull=(TX(TX_idx).Delay/freq); % Angle delays in seconds
                angleDelay=angleDelayFull(tx_apod); % Angle delays of aperture
                forwardDelayAp=Delay(tx_apod); % Forward delays of aperture (no angle)
                [fDelay,fidx]=min(forwardDelayAp); % idx and delay of forward channel
                forwardDelay=fDelay+angleDelay(fidx); % Forward delay (including angle) of the forward piezo
                backwardDelay=Delay(rcv_apod); % Backward delays of just the piezos receiving
                rcvAngles=elementAngles(rcv_apod); % Angles from the pixel to the recieving peizos
                for ch=1:length(backwardDelay) % For each receiving channel
                    timeDelay=(forwardDelay+backwardDelay(ch))'; % Array of all delays from every transmitting piezo to pixel to the receiving piezo
                    sampleDelay=Receive(idx).startSample+floor(timeDelay*freq*2); % Delay in wavelength from the start of the step
                    deltaDelay=2*pi*timeDelay*freq;
                    elementSens=Trans.ElementSens(round((rcvAngles(ch)+pi/2)/(pi/100)+1)); % Sensitivity of recieving piezo for the pixel
                    if elementSens>=Senscutoff % If the sensitivity if over the cutoff
                        if all(sampleDelay>=Receive(idx).startSample) && all(sampleDelay+1<=Receive(idx).endSample) % If its in the sample range
                            RF0=single(RData(sampleDelay,ch)); % Array of all RF0s for every transmitting piezo to pixel to the receiving piezo
                            RF1=single(RData(sampleDelay+1,ch)); % Array of all RF1s for every transmitting piezo to pixel to the receiving piezo
                            cosDelta=cos(deltaDelay);
                            sinDelta=sin(deltaDelay);
                            OutImageIQTmpX=(RF0.*cosDelta+RF1.*sinDelta); % Apply rotation matrix
                            OutImageIQTmpY=(-RF0.*sinDelta+RF1.*cosDelta);
                            if ch > 32*2 && ch <= 32*6 % Flip location of pixel for center 4 rows (not sure why but it fixes a problem)
                                bf(length(X)-x_pix+1,y_pix,z_pix)=bf(length(X)-x_pix+1,y_pix,z_pix)+sum(OutImageIQTmpX)+sqrt(-1)*sum(OutImageIQTmpY); % Add all the data
                            else
                                bf(x_pix,y_pix,z_pix)=bf(x_pix,y_pix,z_pix)+sum(OutImageIQTmpX)+sqrt(-1)*sum(OutImageIQTmpY); % Add all the data
                            end
                        else
                            disp('outside sample range'); % This line of code should never be reached
                            return
                        end
                    end
                end
            end
        end
    end
end
bf=permute(bf,[2 1 3]); % To make same dimensions as Verasonics beamformer