#include "vp_voice_process.h"
#include "vp_aec.h"
#include <math.h>
#include "vp_math.h"

VoiceProcessor::~VoiceProcessor()
{
}

int32_t VoiceProcessor::Initialize(int32_t* param)
{
	cfg_c2r = kiss_fftr_alloc(FftSize_, 0, 0 , 0);
	rx_in = new kiss_fft_scalar[FftSize_];
	cx_out = new kiss_fft_cpx[FftBins_];

	cfg_c2r = kiss_fftr_alloc(FftSize_, 1, 0, 0);
	rx_out = new kiss_fft_scalar[FftSize_];
	cx_in = new kiss_fft_cpx[FftBins_];

	FftWin = new float[FftSize_];
	for (int m = 0; m < FftSize_; ++m) {
		FftWin[m] = sqrtf(0.5f * (1.0f - cosf(2.0f * M_PI * m / (FftSize_ - 1.0f))) );
	}

	FftWinInv = new float[FftSize_];
	for (int m = 0; m < FftSize_; ++m) {
		FftWinInv[m] = 1.0f / (FftWin[m] + 0.02f);
	}

	WinOutComp = new float[FrameSize_];
	int32_t NumAccum = 2 * FftSize_ / FrameSize_;
	float* DataTmp = new float[NumAccum * FrameSize_ + FftSize_];
	VectorClear(DataTmp, NumAccum * FrameSize_ + FftSize_);
	for (int frame = 0; frame <= NumAccum; ++frame) {
		float* DataTmpPoint = DataTmp + frame * FrameSize_;
		for (int m = 0; m < FftSize_; ++m) {
			DataTmpPoint[m] += FftWin[m] * FftWin[m];
		}
	}

	int32_t NumAccum_1 = FftSize_ / FrameSize_ + 1;
	for (int32_t m = 0; m < FrameSize_; ++m) {
		WinOutComp[m] = 1.0f / DataTmp[m + FrameSize_ * NumAccum_1];
	}
	delete[] DataTmp;

	IndexCutDC = Max(1, (int)(150.0f * FftSize_ / SampleRate_ + 0.5f));
	CutDC = new float[IndexCutDC];
	for (int m = 0; m < IndexCutDC; ++m) {
		CutDC[m] = 0.5f * (1.0f - cosf(M_PI * m / (IndexCutDC - 1.0f)));
		CutDC[m] *= CutDC[m];
	}

	// overlap-add
	MicDataAccum = new float[FftSize_];
	VectorClear(MicDataAccum, FftSize_);

	RefDataAccum = new float[FftSize_];
	VectorClear(RefDataAccum, FftSize_);

	OutDataAccum = new float[FftSize_];
	VectorClear(OutDataAccum, FftSize_);

	// frequence data
	MicData = new VP_CPX[FftBins_];
	VectorClear(MicData, FftBins_);

	RefData = new VP_CPX[FftBins_];
	VectorClear(RefData, FftBins_);

	OutData = new VP_CPX[FftBins_];
	VectorClear(OutData, FftBins_);

	vpaec_ = new VPAEC(SampleRate_, FftSize_, FrameSize_);

	Results_ = new InterParamters;
	Results_->ProbFrame = 0.f;
	Results_->EchoStateFlag = false;
	Results_->EchoStateFlagIm = false;
	Results_->ZeroFlag = true;

	Results_->GainEC = new float[FftBins_];
	VectorClear(Results_->GainEC, FftBins_);

	Results_->GainNS = new float[FftBins_];
	VectorClear(Results_->GainNS, FftBins_);

	Results_->GainAGC = 1.f;

	Results_->PsigOut = new float[FftBins_];
	VectorClear(Results_->PsigOut, FftBins_);

	Results_->PsigIn = new float[FftBins_];
	VectorClear(Results_->PsigIn, FftBins_);

	Results_->Psig = new float[FftBins_];
	VectorClear(Results_->Psig, FftBins_);

	Results_->Pnoise = Results_->PsigOut;

}

void VoiceProcessor::Release()
{
	free(cfg_r2c);
	delete[] rx_in;
	delete[] cx_out;

	free(cfg_c2r);
	delete[] rx_out;
	delete[] cx_in;

	delete[] FftWin;
	delete[] WinOutComp;
	delete[] FftWinInv;
	
	delete[] CutDC;

	delete[] MicData;
	delete[] RefData;
	delete[] OutData;

	delete[] MicDataAccum;
	delete[] RefDataAccum;
	delete[] OutDataAccum;

	vpaec_->Release();
	delete vpaec_;

	delete[] Results_->GainEC;
	delete[] Results_->GainNS;
	delete[] Results_->Psig;
	delete[] Results_->PsigIn;
	delete[] Results_->PsigOut;
	delete Results_;

}

void VoiceProcessor::Process(int16_t* micdata, int16_t* refdata, int16_t* outdata)
{
	if (micdata == NULL || refdata == NULL || outdata == NULL) {
		return;
	}

	int32_t num = FrameSize_;

	int NumFrame = num / FrameSize_;
	for (int32_t kframe = 0; kframe < NumFrame; ++kframe) {
		int16_t* outdata_point = outdata + kframe * FrameSize_;
		int16_t* micdata_point = micdata + kframe * FrameSize_;
		int16_t* refdata_point = refdata + kframe * FrameSize_;

		float* MicDataAccumPoint = MicDataAccum + FftSize_ - FrameSize_;
		float* RefDataAccumpoint = RefDataAccum + FftSize_ - FrameSize_;
		for (int m = 0; m < FrameSize_; ++m) {
			MicDataAccumPoint[m] = (float)micdata_point[m];
			RefDataAccumpoint[m] = (float)refdata_point[m];
		}

		/* Fft for mic data */
		for (int32_t m = 0; m < FftSize_; ++m) {
			rx_in[m] = MicDataAccum[m] * FftWin[m];
		}
		kiss_fftr(cfg_r2c, rx_in, cx_out);
		for (int32_t m = 0; m < FftBins_; ++m) {
			MicData[m] = VP_CPX(cx_out[m].r, cx_out[m].i);
		}

		/* Fft for ref data */
		for (int32_t m = 0; m < FftSize_; ++m) {
			rx_in[m] = RefDataAccum[m] * FftWin[m];
		}
		kiss_fftr(cfg_c2r, rx_in, cx_out);
		for (int32_t m = 0; m < FftBins_; ++m) {
			RefData[m] = VP_CPX(cx_out[m].r, cx_out[m].i);
		}

		FrameProcess();

		/* IFFT */
		for (int32_t m = 0; m < IndexCutDC; ++m) {
			cx_in[m].r = OutData[m].re() * CutDC[m];
			cx_in[m].i = OutData[m].im() * CutDC[m]; 
		}
		for (int32_t m = IndexCutDC; m < FftBins_; ++m) {
			cx_in[m].r = OutData[m].re();
			cx_in[m].i = OutData[m].im();
		}
		kiss_fftri(cfg_c2r, cx_in, rx_out);

		for (int32_t m = 0; m < FftSize_; ++m) {
			OutDataAccum[m] += rx_out[m] * FftWin[m];
		}

		for (int32_t m = 0; m < FrameSize_; ++m) {
			outdata_point[m] = (int16_t)(OutDataAccum[m] * WinOutComp[m]);
		}
		
		memmove(OutDataAccum, OutDataAccum + FrameSize_, 4 * (FftSize_ - FrameSize_));
		memmove(MicDataAccum, MicDataAccum + FrameSize_, 4 * (FftSize_ - FrameSize_));
		memmove(RefDataAccum, RefDataAccum + FrameSize_, 4 * (FftSize_ - FrameSize_));

		memset(OutDataAccum + FftSize_ - FrameSize_, 0, FrameSize_ * 4);
	}
}

void VoiceProcessor::FrameProcess()
{
	vpaec_->Process(Results_, MicData, RefData, OutData);

	// vpans
}