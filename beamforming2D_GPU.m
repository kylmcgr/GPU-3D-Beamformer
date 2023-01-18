
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                             %
%     Power Doppler acquisition script        %
%                                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear, format compact

% set global variabls & paths
global UF na path_save dir_save 
dir_save = ['D:\Data_GPU_fUS\Test\'];
mkdir(dir_save)

% UF acquisition parameters & Image parameters
UF.Probe            = 'L22-14vX';
UF.Depth(1)         = 2.5;           % initial depth [mm]
UF.Depth(2)         = 9.5;          % final depth [mm]
UF.NbOfBlocs        = 2;           % number of loops
UF.TwFreq           = 15.625;      % emission frequency
UF.RcvFreq          = UF.TwFreq;   % receive frequency
UF.DutyCycle        = 0.67;        % 0.75 or 1 for HF
UF.NbHalfCycle      = 3;           % Nb of Half cycle in emission
UF.Time_loop        = 1;           % [s]
UF.ImgVoltage       = 25;          % (V) Voltage applied to the transducers     


% UF block parameters

UF.numFrames        = 220;
UF.dopAngle         = (-14:14/7:14) * pi/180;
na                  = length(UF.dopAngle);
UF.dopFrameRate     = 500;
UF.sampling_mode    = 2;  % 2 = 100%, 4 = 200%

UF.AntiAliasingFilter  = 20;

%% save UF params 
UF.path_save = path_save;
UF.dir_save = dir_save;

save([dir_save 'UF'],'UF')

%% Define system parameters.
Resource.Parameters.numTransmit = 128;      % number of transmit channels.
Resource.Parameters.numRcvChannels = 128;    % number of receive channels.
Resource.Parameters.speedOfSound = 1540;    % set speed of sound in m/sec before calling computeTrans
Resource.Parameters.verbose = 2;
Resource.Parameters.initializeOnly = 0;
Resource.Parameters.simulateMode = 0;
Resource.Parameters.Connector = 1;

Trans.name = UF.Probe;
Trans.units = 'wavelengths';    % Explicit declaration avoids warning message when selected by default
Trans.frequency = UF.TwFreq ;   % emission frequency [MHz]
Trans = computeTrans(Trans);  % LA-16 transducer is 'known'(added by user) transducer so we can use computeTrans.
Trans.maxHighVoltage = 25;      % set maximum high voltage limit for pulser supply.

% conversion in lambda
UF.Lambda = Resource.Parameters.speedOfSound/UF.RcvFreq*1e-3;   % [mm]
P.startDepth =  floor(UF.Depth(1)/UF.Lambda);	% Acquisition depth in wavelengths
P.endDepth =    ceil(UF.Depth(2)/UF.Lambda);    % This should preferrably be a multiple of 128 samples

% Specify PData structure array.
PData(1).PDelta = [Trans.spacing, 0, 1];
PData(1).Size(1) = ceil((P.endDepth-P.startDepth)/PData(1).PDelta(3));
PData(1).Size(2) = ceil((Trans.numelements*Trans.spacing)/PData(1).PDelta(1));
PData(1).Size(3) = 1;      % single image page
PData(1).Origin = [-Trans.spacing*(Trans.numelements-1)/2,0,P.startDepth]; % x,y,z of upper lft crnr.

% Specify Resources.
maxAcqLength = ceil(sqrt(P(1).endDepth^2 + ((Trans.numelements-1)*Trans.spacing)^2)); % this will be used in the receive struct
Resource.RcvBuffer(1).datatype = 'int16';
Resource.RcvBuffer(1).rowsPerFrame = 128*ceil(maxAcqLength/128*2)*na*UF.sampling_mode; % this size allows for maximum range %modif
Resource.RcvBuffer(1).colsPerFrame = Resource.Parameters.numRcvChannels;
Resource.RcvBuffer(1).numFrames = UF.numFrames;    % Nb of frames stored in RcvBuffer.

% Interbuffer used to store IQ data
Resource.InterBuffer(1).datatype = 'complex';
Resource.InterBuffer(1).numFrames = 1;                      % one intermediate buffer needed.
Resource.InterBuffer(1).pagesPerFrame = na*UF.numFrames;
Resource.InterBuffer(1).rowsPerFrame = PData(1).Size(1);	% this size allows for maximum range
Resource.InterBuffer(1).colsPerFrame = PData(1).Size(2);

Resource.ImageBuffer(1).datatype = 'double';
Resource.ImageBuffer(1).numFrames = 1;
Resource.ImageBuffer(1).rowsPerFrame = PData(1).Size(1);
Resource.ImageBuffer(1).colsPerFrame = PData(1).Size(2);

% Specify Transmit waveform structure.  
TW(1).type = 'parametric';
TW(1).Parameters = [UF.TwFreq,UF.DutyCycle,UF.NbHalfCycle,1];

UF.TW = TW;

% Specify TX structure array.
TX = repmat(struct('waveform', 1, ...
                   'Origin', [0.0,0.0,0.0], ...
                   'Apod', kaiser(Resource.Parameters.numTransmit,1)', ...
                   'focus', 0.0, ...
                   'Steer', [0.0,0.0], ...
                   'Delay', zeros(1,Trans.numelements)), 1, na);
               
% - Set event TX attributes.
for n = 1:na   % na transmit events
    TX(n).Steer(1) = UF.dopAngle(n);
    TX(n).Delay = computeTXDelays(TX(n));
end

UF.TX = TX;

% Specify TGC Waveform structure.
TGC.CntrlPts =  512*ones(1,8); % [330 560 780 1010 1023 1023 1023 1023];
TGC.rangeMax = P.endDepth;
TGC.Waveform = computeTGCWaveform(TGC);

% Specify Receive structure arrays.  

% For Doppler, use narrow bandwidth coefficients (50% BW) centered at
% 15.625; this is a copy of update function's default coef array for 1
% samples per wave
BPFDop = [ -0.00162 +0.00000 +0.00568 +0.00000 -0.01065 +0.00000 +0.01349 ...
           +0.00000 -0.00858 +0.00000 -0.00955 +0.00000 +0.04312 +0.00000 ...
           -0.08841 +0.00000 +0.13550 +0.00000 -0.17130 +0.00000 +0.18463 ];

switch UF.sampling_mode
    case 4
        sampleMode_str = 'NS200BW';
        Receive = repmat(struct('Apod', ones(1,Trans.numelements), ...
                        'startDepth', P.startDepth, ...
                        'endDepth', maxAcqLength,...
                        'bufnum', 1, ...
                        'framenum', 1, ...
                        'acqNum', 1, ...
                        'TGC', 1, ...
                        'InputFilter', BPFDop,...
                        'sampleMode', sampleMode_str, ...
                        'mode', 0), ...
                        1, na*Resource.RcvBuffer(1).numFrames);
    case 2
        sampleMode_str = 'BS100BW';
        Receive = repmat(struct('Apod', ones(1,Trans.numelements), ...
                        'startDepth', P.startDepth, ...
                        'endDepth', maxAcqLength,...
                        'bufnum', 1, ...
                        'framenum', 1, ...
                        'acqNum', 1, ...
                        'TGC', 1, ...
                        'InputFilter', BPFDop,...
                        'sampleMode', sampleMode_str, ...
                        'demodFrequency', Trans.frequency, ...
                        'mode', 0), ...
                        1, na*Resource.RcvBuffer(1).numFrames);
end
          

% - Set event specific Receive attributes for each frame.
for i = 1:Resource.RcvBuffer(1).numFrames
    for j = 1:na
        
        Receive(na*(i-1)+j).framenum = i;
        Receive(na*(i-1)+j).acqNum = j;
    end
end

%% Define ReconInfo structures.
% We need na ReconInfo structures for na steering angles.
ReconInfo = repmat(struct('mode', 'replaceIQ', ...  % default is to accumulate IQ data.
                   'txnum', 1, ...
                   'rcvnum', 1, ...
                   'pagenum',[],...
                   'regionnum', 1), 1, na);

%_________________________________________________________________________%
% Recon/ReconInfo for reconstruct and store the IQData of each acq        %

% Specify Recon structure arrays.
% - We need one Recon structures which will be used for each frame.
k = 0;
page = 1 ;
for frame = 1:UF.numFrames
    Recon(frame) = struct('senscutoff', 0.8, ...
        'pdatanum', 1, ...
        'rcvBufFrame',frame, ...
        'IntBufDest', [1,1], ...
        'ImgBufDest', [0,0], ...
        'newFrameTimeout', UF.Time_loop*1e3,...  % needed for the software to wait for a new frame in the first rcv slot
        'RINums',(k+1:k+na));
    k = k+na;
    
    % Define ReconInfo structures.
    % We need na ReconInfo structures for na steering angles.
    
    for angl = 1:na
        ReconInfo(na*(frame-1) + angl) = struct('mode', 'replaceIQ', ...   % replace IQ data every time
            'txnum', angl, ...
            'rcvnum', na*(frame-1) + angl, ...       % the index of the receive acquisition
            'pagenum', page,...
            'regionnum', 1);
        page = page + 1 ;
    end
end


%_________________________________________________________________________%
%                                                                         %

% Specify Process structure array.
% Filter IQ data

Process(1).classname = 'External';     % process structure for 1st Doppler ensemble
Process(1).method = 'computeSVD';
Process(1).Parameters = {'srcbuffer','inter',... % name of buffer to process.
                         'srcbufnum',1,...
                         'srcframenum',1,...            % process the most recent frame.
                         'dstbuffer','image',...
                         'dstbufnum', 1,...
                         'dstframenum',1};
EF(1).Function = vsv.seq.function.ExFunctionDef('computeSVD', @computeSVD);

                     
% Save IQ data
Process(2).classname = 'External';     % process structure for 1st Doppler ensemble
Process(2).method = 'saveIQData';
Process(2).Parameters = {'srcbuffer','inter',... % name of buffer to process.
                         'srcbufnum',1,...
                         'srcframenum',1,...            % process the most recent frame.
                         'dstbuffer','none'};
EF(2).Function = vsv.seq.function.ExFunctionDef('saveIQData', @saveIQData);


% Display block time
Process(3).classname = 'External';     % process structure for 1st Doppler ensemble
Process(3).method = 'dispBlockTime';
Process(3).Parameters = {'srcbuffer','none',... 
                         'dstbuffer','none'};
EF(3).Function = vsv.seq.function.ExFunctionDef('dispBlockTime', @dispBlockTime);


%% Specify SeqControl structure arrays. 
maxReturnTrip = sqrt(P(1).endDepth^2 + (Trans.numelements*Trans.spacing)^2)*2;
maxTOF = round(2*maxReturnTrip/Trans.frequency);  % to be secure, 2 times
time_compound = na*maxTOF;
time_ensemble = 1/(UF.dopFrameRate*1e-6)-(time_compound-maxTOF);
time_2_nextSEQ = UF.Time_loop*1e6- 1/(UF.dopFrameRate*1e-6)*UF.numFrames - time_ensemble;

% Specify SeqControl structure arrays.
SeqControl(1).command = 'jump'; % jump back to start, not used here
SeqControl(1).argument = 1;
SeqControl(2).command = 'timeToNextAcq';  % time between synthetic aperture acquisitions
SeqControl(2).argument = maxTOF;  % 200 usecs
SeqControl(3).command = 'triggerOut';
SeqControl(4).command = 'returnToMatlab';
SeqControl(5).command = 'timeToNextAcq';  % time between frames
SeqControl(5).argument = time_ensemble;  % 2 msec (<=> 500 Hz)
SeqControl(6).command = 'timeToNextAcq';  % time between blocks
SeqControl(6).argument = time_2_nextSEQ;  % 3 sec 
SeqControl(7).command = 'sync';  % synchronisation soft hard, (no ping pong)
SeqControl(7).argument = UF.Time_loop*1e6;  % time out for sync  % needed for the software to wait for the hardware to end the last TTNA of the loop
SeqControl(8).command = 'loopCnt'; % - Set loop count. for looping on blocs
SeqControl(8).argument = UF.NbOfBlocs-1;  %
SeqControl(8).condition = 'counter1';
SeqControl(9).command = 'loopTst';  % loop test
SeqControl(9).argument = [];    % set apres
SeqControl(9).condition = 'counter1';

nsc = 10; % nsc is count of SeqControl objects

% Specify Event structure arrays.
n = 1;

% set loop count
Event(n).info = 'start counter';
Event(n).tx = 0;   % use next TX structure.
Event(n).rcv = 0;
Event(n).recon = 0;      % no reconstruction.
Event(n).process = 0;    % no processing
Event(n).seqControl = 8;
n = n+1;
SeqControl(9).argument = n;

for i = 1:Resource.RcvBuffer(1).numFrames
    for j = 1:na                    % Acquire frame
            Event(n).info = 'Acquire RF';
            Event(n).tx = j;   % use next TX structure.
            Event(n).rcv = na*(i-1)+j ;
            Event(n).recon = 0;      % no reconstruction.
            Event(n).process = 0;    % no processing
            Event(n).seqControl = 2;
            n = n+1;
    end
    
    % set sync hardware and software for the first TX
    Event(2).seqControl = [7,2];
    
    % Set last acquisitions SeqControl for transferToHost.
    Event(n-1).seqControl = [5,nsc]; % modify last acquisition Event's seqControl
      SeqControl(nsc).command = 'transferToHost'; % transfer all acqs to host buffer
      nsc = nsc+1;

    Event(n).info = 'Reconstruct'; 
    Event(n).tx = 0;        % no transmit
    Event(n).rcv = 0;       % no rcv
    Event(n).recon = i;     % reconstruction
    Event(n).process = 0;	% process
    Event(n).seqControl = 0;
    n = n+1;
end

Event(n-2).seqControl = [6,nsc-1]; % modify last acquisition of the lat frame TTNA

Event(n).info = 'Sequence Time Control'; 
Event(n).tx = 0;         % no transmit
Event(n).rcv = 0;        % no rcv
Event(n).recon = 0;      % no reconstruction
Event(n).process = 3;    % process
Event(n).seqControl = 0;
n = n+1;

Event(n).info = 'save IQ Data Process'; 
Event(n).tx = 0;         % no transmit
Event(n).rcv = 0;        % no rcv
Event(n).recon = 0;      % no reconstruction
Event(n).process = 2;    % process
Event(n).seqControl = 0;
n = n+1;

Event(n).info = 'loop everything'; 
Event(n).tx = 0;         % no transmit
Event(n).rcv = 0;        % no rcv
Event(n).recon = 0;      % no reconstruction
Event(n).process = 0;    % process
Event(n).seqControl = 9 ;
n = n+1;

% Rcv profile
RcvProfile.DCsubtract = 'on';               % substract DC signal if 'on'
RcvProfile.AntiAliasCutoff = UF.AntiAliasingFilter;    % antialiasing analogical filter cuttoff freq [MHz]
RcvProfile.PgaGain = 30;	% 24 ou 30	% analog gain in dBfor preamp  %#test
RcvProfile.LnaGain = 24;	% gain in dB of the fixed gain low noise amp
RcvProfile.LnaZinSel = 31;	% Force high-Z state for best Doppler sensitivity  %#test



%Check coherence of the VS structures

if max([Event(:).rcv])>length(Receive)
    warning('Probably a problem of receive indexing in the Event struct')
elseif max([Event(:).tx])>length(TX)
    warning('Probably a problem of TX indexing in the Event struct')
elseif max([Event(:).seqControl])>length(SeqControl)
    warning('Probably a problem of SeqControl indexing in the Event struct')
elseif max([ReconInfo(:).txnum])>length(TX)
    warning('Probably a problem of TX indexing in the ReconInfo struct')
% elseif max([ReconInfo(:).pagenum])>UF.numFrames
%     warning('Probably a problem of interbuffer pages management in the ReconInfo struct')
elseif max([ReconInfo(:).rcvnum])>length(Receive)
    warning('Probably a problem of receive indexing in the ReconInfo struct')
end

% UI controls

% Set TPCHighVoltage for profile one to 20V
UI(1).Statement = '[result,hv] = setTpcProfileHighVoltage(UF.ImgVoltage,1);';
UI(2).Statement = 'hv1Sldr = findobj(''Tag'',''hv1Sldr'');';
UI(3).Statement = 'set(hv1Sldr,''Value'',hv);';
UI(4).Statement = 'hv1Value = findobj(''Tag'',''hv1Value'');';
UI(5).Statement = 'set(hv1Value,''String'',num2str(hv,''%.1f''));';


% Save all the structures to a .mat file. & auto start
filename = 'MatFiles\fUSloop';
save(filename);
save([dir_save 'UF'],'UF');


VSX
return
% 
%% External function definitions

function Dop = computeSVD(RData)
    global na;
    IQ = squeeze(RData);
    IQ = reshape(IQ, size(IQ, 1), size(IQ, 2), na, []);
    IQ = squeeze(mean(IQ, 3));
  
    ncut = 30;
    % SVD PowerDoppler
    IQ_signal = IQ;
    [nz, nx, nt] = size(IQ_signal);
    IQ_signal = reshape(IQ_signal, [nz*nx, nt]);
    cov_matrix = IQ_signal'*IQ_signal;
    [Eig_vect, Eig_val]= eig(cov_matrix);
    Eig_vect=fliplr(Eig_vect);
    Eig_val=rot90(Eig_val,2);
    M_A = IQ_signal*Eig_vect;
    skipped_eig_val =[1:ncut];
    IQF_tissu = M_A(:,skipped_eig_val)*Eig_vect(:,skipped_eig_val)';
    IQF_tissu = reshape(IQF_tissu, [nz, nx, nt]);
    IQ_signal = reshape(IQ_signal, [nz, nx, nt]);
    IQF_corrected = IQ_signal-IQF_tissu;
    
    Dop = mean(abs(IQF_corrected(:,:,:)).^2,3);
    Dop = (Dop - min(Dop(:)))./range(Dop(:));
    
    figure(56)
    imagesc(interp2(sqrt(Dop),2)), colormap hot(256)
    caxis([0.05 0.95])
end



function saveIQData(IData,QData)
IQData=complex(IData,QData);

    persistent bloc_count;
    global na UF dir_save;
    if isempty(bloc_count)
        bloc_count = 1;
    end
    tmp = squeeze(IQData);
    tmp = reshape(tmp, size(tmp, 1), size(tmp, 2), na, []);
    tmp = squeeze(mean(tmp, 3));   % do the compounding (less disk space needed)    
    IQ = [real(tmp) imag(tmp)];
    file_name = sprintf('fUS_block_%.3d.bin', bloc_count);
    fid = fopen([dir_save file_name],'W');   % fast save
    fwrite(fid,IQ, 'double');
    fclose(fid);
    bloc_count = bloc_count+1;
end




function dispBlockTime()
    persistent time_bloc;
    global dir_save;
    if isempty(time_bloc)
        time_bloc = tic;
    end
    T = toc(time_bloc) % absolute time (ms)
    D = datestr(now, 'HH-MM-SS.FFF'); % fUS image interval (s)
    fidresult = fopen([dir_save '\timeLog.txt'], 'a');
    fprintf(fidresult,'%s  %g\r\n', D, T);
    fclose(fidresult);
end
