#include "kiss_fft.h"
#include "kiss_fftr.h"
#include "_kiss_fft_guts.h"
#include "vp_math.h"
#include "vp_cpx.h"
#include "vp_aec.h"

class VPAEC;

class VoiceProcessor
{
private:
	const int32_t SampleRate_;
	const int32_t FftSize_;
	const int32_t FftBins_;
	const int32_t FrameSize_;

public:
	VoiceProcessor(int32_t sampelrate, int32_t fftsize, int32_t framesize)
		: SampleRate_(sampelrate), FftSize_(fftsize), FftBins_(fftsize / 2 + 1), FrameSize_(framesize) {}

	~VoiceProcessor();

	int32_t Initialize(int32_t* params);

	void Release();

	void Process(int16_t* micdata, int16_t* refdata, int16_t* outdata);

private:
	void FrameProcess();

	kiss_fftr_cfg cfg_r2c;
	kiss_fft_cpx *cx_out;
	kiss_fft_scalar *rx_in;

	kiss_fftr_cfg cfg_c2r;
	kiss_fft_cpx *cx_in;
	kiss_fft_scalar *rx_out;

	float* FftWin;
	float* WinOutComp;
	float* FftWinInv;

	float* CutDC;
	int32_t IndexCutDC;

	float* MicDataAccum;
	float* RefDataAccum;
	float* OutDataAccum;

	VP_CPX* MicData;
	VP_CPX* RefData;
	VP_CPX* OutData;

	VPAEC* vpaec_;

	InterParamters* Results_;

};