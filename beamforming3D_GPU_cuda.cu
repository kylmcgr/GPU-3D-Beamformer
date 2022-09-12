// GPU beamformer CUDA code for 3D IQ data
// Beamforms a single pixel, called by the beamforming kernel
// Author - Kyle McGraw
// Written on August 2022
// 

__constant__ float invSpeedOfSound;
__constant__ int xPixLen;
__constant__ int yPixLen;
__constant__ int zPixLen;
__constant__ int xPiezoLen;
__constant__ int yPiezoLen;
__constant__ int receiveLen;
__constant__ float freq;
__constant__ float pi;
__constant__ int sampleDim;
__constant__ float sensCutoff;
__constant__ float elementSens[101];

__global__ void beamforming3D_GPU_cuda(
	float2* bf,
    const short* __restrict__ RData,
    const float* __restrict__ grid_X,
    const float* __restrict__ grid_Y,
    const float* __restrict__ grid_Z,
    const float* __restrict__ piezo_X,
    const float* __restrict__ piezo_Y,
    const int* __restrict__ idxs_TX,
    const int* __restrict__ TX_apertures,
    const int* __restrict__ RCV_apertures,
    const float* __restrict__ angleDelays,
    const int* __restrict__ startSamples,
    const int* __restrict__ endSamples) {
    const unsigned int x_pix = blockIdx.x * blockDim.x + threadIdx.x; // x coord of the pixel
    const unsigned int y_pix = blockIdx.y * blockDim.y + threadIdx.y; // y coord of the pixel
    const unsigned int z_pix = blockIdx.z * blockDim.z + threadIdx.z; // z coord of the pixel
    if  (x_pix < xPixLen && y_pix < yPixLen && z_pix < zPixLen) {  // make sure its a real pixel because round up the number of thread blocks
        for (int idx = 0; idx < receiveLen; idx++) { // For each transmit-receive step
            const unsigned int TX_idx = idxs_TX[idx]; // Index in TX object
            const unsigned int tx_aperture = TX_apertures[TX_idx]; // Transmitting aperture
            const unsigned int rcv_aperture = RCV_apertures[idx]; // Receiving aperture
            for (int xr = 0; xr < xPiezoLen; xr++) { // For each receiving piezo
                for (int yr = (rcv_aperture-1)*8; yr < (rcv_aperture-1)*8+8; yr++) {
                    const float returnDelay = sqrtf((grid_X[x_pix]-piezo_X[xr])*(grid_X[x_pix]-piezo_X[xr])+(grid_Y[y_pix]-piezo_Y[yr])*(grid_Y[y_pix]-piezo_Y[yr])+grid_Z[z_pix]*grid_Z[z_pix])*invSpeedOfSound;
                    const unsigned int channelr = (yr-((rcv_aperture-1)*8))*32+xr; // Receiving channel 0-255, used for indexing RData channels
                    const unsigned int sensIdx = round((atanf(sqrtf((grid_X[x_pix]-piezo_X[xr])*(grid_X[x_pix]-piezo_X[xr])+(grid_Y[y_pix]-piezo_Y[yr])*(grid_Y[y_pix]-piezo_Y[yr]))/grid_Z[z_pix])+pi/2)/(pi/100));
                    const float piezoSens=elementSens[sensIdx];
                    if (piezoSens>=sensCutoff) { // If the sensitivity if over the cutoff
                        int minXt = 0; // Find the closest tranmitting piezo
                        int minYt = (tx_aperture-1)*8;
                        for (int xt = 0; xt < xPiezoLen; xt++) { // For each transmitting piezo
                            if (fabsf(grid_X[x_pix]-piezo_X[xt]) < fabsf(grid_X[x_pix]-piezo_X[minXt])){
                                minXt=xt;
                            }
                        }
                        for (int yt = (tx_aperture-1)*8; yt < (tx_aperture-1)*8+8; yt++) {
                            if (fabsf(grid_Y[y_pix]-piezo_Y[yt]) < fabsf(grid_Y[y_pix]-piezo_Y[minYt])){
                                minYt=yt;
                            }
                        }
                        const float forwardDelay = sqrtf((grid_X[x_pix]-piezo_X[minXt])*(grid_X[x_pix]-piezo_X[minXt])+(grid_Y[y_pix]-piezo_Y[minYt])*(grid_Y[y_pix]-piezo_Y[minYt])+grid_Z[z_pix]*grid_Z[z_pix])*invSpeedOfSound;;
                        const unsigned int channelt = minYt*32+minXt; // Transmitting channel 0-1023, used for indexing angle delay
                        const float angleDelay = angleDelays[TX_idx*1024+channelt]/freq; // Angle delay in seconds from linearlized 2D array
                        const float timeDelay = forwardDelay+returnDelay+angleDelay; // Delay from transmitting piezo to pixel to receiving piezo
                        const unsigned int sampleDelay = startSamples[idx]+floor(timeDelay*freq*2); // Delay in wavelength from the start of the step
                        const float deltaDelay=2*pi*timeDelay*freq;
                        if (sampleDelay>=startSamples[idx] && sampleDelay+1<=endSamples[idx]){ // If its in the sample range
                            const short RF0 = RData[sampleDelay+sampleDim*channelr]; // Get the data from linearlized 2D array
                            const short RF1 = RData[sampleDelay+1+sampleDim*channelr];
                            const float cosDelta = cosf(deltaDelay);
                            const float sinDelta = sinf(deltaDelay);
                            if (yr >= (rcv_aperture-1)*8+2 && yr < (rcv_aperture-1)*8+6){
                                bf[xPixLen-1-x_pix+y_pix*xPixLen+z_pix*xPixLen*yPixLen].x += RF0*cosDelta+RF1*sinDelta; // Apply rotation matrix
                                bf[xPixLen-1-x_pix+y_pix*xPixLen+z_pix*xPixLen*yPixLen].y += -RF0*sinDelta+RF1*cosDelta;
                            } else {
                                bf[x_pix+y_pix*xPixLen+z_pix*xPixLen*yPixLen].x += RF0*cosDelta+RF1*sinDelta; // Apply rotation matrix
                                bf[x_pix+y_pix*xPixLen+z_pix*xPixLen*yPixLen].y += -RF0*sinDelta+RF1*cosDelta;
                            }
                        }
                    }
                }
            }
        }
    }
}