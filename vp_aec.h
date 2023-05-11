#include "vp_cpx.h"
#include "vp_typedef.h"

class VPAEC
{
public:
	const int32_t SampleRate_;
	const int32_t FftSize_;
	const int32_t FftBins_;
	const int32_t FrameSize_;

public:
	VPAEC(int32_t samplerate, uint32_t fftsize, uint32_t framesize)
		:SampleRate_(samplerate), FftSize_(fftsize), FrameSize_(framesize), FftBins_(fftsize / 2 + 1) {}

	~VPAEC();

	int32_t Initialize(int32_t* param);

	void Release();

	void Process(InterParamters* Params, const VP_CPX* MicData, const VP_CPX* RefData, VP_CPX* OutData);

private:
	int32_t AECValid;
	int32_t EC_REF_ALIGN_MAX;
	
	int32_t MICNUM;
	int32_t REFNUM;
	int32_t SEARCHNUM;
	int32_t SMOOTHBINS;
	
	float* WinSmooth;

	VP_CPX* MicDataBuffer;
	VP_CPX* RefDataBuffer;
	float* MicDataBufferAbs2;
	float* RefDataBufferAbs2;

	float* MicDataAbs2;
	float* RefDataAbs2;

	int32_t IndexMic;
	int32_t IndexRef;

	int32_t GetMicIndex(int32_t i);
	int32_t GetRefIndex(int32_t i);

	Feature* MicFeature;

	Feature* RefFeature;

	float* FeatureCurrent;
	float* FeatureSmooth;
	float* FeatureCopy;

	int32_t BinStart;
	int32_t BinStep;

	const float* GetFeatureCurrent()const;
	const float* GetFeatureSmooth()const;

	int32_t Bcounter(uint32_t n);

	Feature FeatureCaculate(const float* data2, const float bias);
	
	void FeatureCurrCalcu();

	void FeatureSmthCalcu();

	int32_t SmoothCounter;
	int32_t IndexAlignment;

	/* parameters of adaptive filters */
	float* DetectionFlag;
	VP_CPX* Wa;
	VP_CPX* Pmr;
	float* Prr;
	int32_t NumStage_a;

	int32_t* IndexMic4EC;
	int32_t* IndexRef4EC;

	void DetectModule();
	float* CorrTmp1;
	float* CorrTmp2;
	VP_CPX* CorrTmpc;
	int32_t DetectionBins_s;
	int32_t DetectionBins_e;

	float ThreshBand;
	float ThreshBins;

	void EchoCanceling(InterParamters* Params, VP_CPX* OutData);
	int32_t IndexLowFreq;
	int32_t IndexMidFreq;

	void EchoCancelingProcess(VP_CPX* OutData, VP_CPX* EchoEst, int IndexAlign);
	VP_CPX* EchoEst;
	VP_CPX* EchoEst_bk;
	VP_CPX* OutData_bk;

	void EchoCancelingProcess4Search();
	int32_t SearchStart;
	int32_t SearchEnd;
	int32_t SearchNum;
	int32_t NumStage_s;
	VP_CPX* Wa_s;
	VP_CPX* Pmr_s;
	float* Prr_s;
	int32_t* SearchFrameValid;
	float* CancelOutSumCurr;
	float* CancelOutSumSmth;

	void DetectModule4Search();
	int32_t IndexAlignSearch;

	void EchoStateCalculate(InterParamters* Params);
	float Pecho4ES;
	float Psig4ES;
	int32_t Counter4ES;
	int32_t Index4ES_s;
	int32_t Index4ES_e;
	float PechoSmooth;
	float PsigSmooth;
	float SEROut;

	void EchoSuppression(InterParamters* Params);
	float* EchoEstAbs2;
	float* Ps;
	float* Pe;
	float* Ps0;
	float* Pe0;
	float* Gain0;
	int32_t PostBins_s;
	int32_t PostBins_e;
	int32_t PostState;
	int32_t PostStateCounter;

	int32_t EchoStateCounter;
	float* PResidualEcho;
	float* PResidualEcho0;


};












