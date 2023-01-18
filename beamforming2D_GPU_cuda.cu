// GPU beamformer CUDA code for 2D IQ data
// Beamforms a single pixel, called by the beamforming kernel
// Author - Kyle McGraw
// Written on January 2023
// 

__constant__ float invSpeedOfSound;
__constant__ int xPixLen;
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
    const float* __restrict__ grid_Z,
    const float* __restrict__ piezo_X,
    const int* __restrict__ idxs_TX,
    const int* __restrict__ TX_apertures,
    const int* __restrict__ RCV_apertures,
    const float* __restrict__ angleDelays,
    const int* __restrict__ startSamples,
    const int* __restrict__ endSamples) {
    const unsigned int x_pix = blockIdx.x * blockDim.x + threadIdx.x; // x coord of the pixel
    const unsigned int z_pix = blockIdx.y * blockDim.y + threadIdx.y; // z coord of the pixel
    if  (x_pix < xPixLen && z_pix < zPixLen) {  // make sure its a real pixel because round up the number of thread blocks
        for (int idx = 0; idx < receiveLen; idx++) { // For each transmit-receive step
            const unsigned int TX_idx = idxs_TX[idx]; // Index in TX object
            const unsigned int tx_aperture = TX_apertures[TX_idx]; // Transmitting aperture
            const unsigned int rcv_aperture = RCV_apertures[idx]; // Receiving aperture
            for (int xr = 0; xr < xPiezoLen; xr++) { // For each receiving piezo
                const float returnDelay = sqrtf((grid_X[x_pix]-piezo_X[xr])*(grid_X[x_pix]-piezo_X[xr])+grid_Z[z_pix]*grid_Z[z_pix])*invSpeedOfSound;
                const unsigned int sensIdx = round((atanf((grid_X[x_pix]-piezo_X[xr])/grid_Z[z_pix])+pi/2)/(pi/100));
                const float piezoSens=elementSens[sensIdx];
                if (piezoSens>=sensCutoff) { // If the sensitivity if over the cutoff
                    int minXt = 0; // Find the closest tranmitting piezo
                    for (int xt = 1; xt < xPiezoLen; xt++) { // For each transmitting piezo
                        if (fabsf(grid_X[x_pix]-piezo_X[xt]) < fabsf(grid_X[x_pix]-piezo_X[minXt])){
                            minXt=xt;
                        }
                    }
                    const float forwardDelay = sqrtf((grid_X[x_pix]-piezo_X[minXt])*(grid_X[x_pix]-piezo_X[minXt])+grid_Z[z_pix]*grid_Z[z_pix])*invSpeedOfSound;;
                    const float angleDelay = angleDelays[TX_idx*128+minXt]/freq; // Angle delay in seconds from linearlized 2D array
                    const float timeDelay = forwardDelay+returnDelay+angleDelay; // Delay from transmitting piezo to pixel to receiving piezo
                    const unsigned int sampleDelay = startSamples[idx]+floor(timeDelay*freq*2); // Delay in wavelength from the start of the step
                    const float deltaDelay=2*pi*timeDelay*freq;
                    if (sampleDelay>=startSamples[idx] && sampleDelay+1<=endSamples[idx]){ // If its in the sample range
                        const short RF0 = RData[sampleDelay+sampleDim*xr]; // Get the data from linearlized 2D array
                        const short RF1 = RData[sampleDelay+1+sampleDim*xr];
                        const float cosDelta = cosf(deltaDelay);
                        const float sinDelta = sinf(deltaDelay);
                        bf[x_pix+z_pix*xPixLen].x += RF0*cosDelta+RF1*sinDelta; // Apply rotation matrix
                        bf[x_pix+z_pix*xPixLen].y += -RF0*sinDelta+RF1*cosDelta;
                    }
                }
            }
        }
    }
}