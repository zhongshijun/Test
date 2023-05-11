#include "vp_cpx.h"
#include "vp_math.h"
#include "vp_typedef.h"

class VPANS
{
private:
	const int32_t SampleRate_;
	const int32_t FftSize_;
	const int32_t FftBins_;
	const int32_t FrameSize_;

public:
	VPANS(int32_t samplerate, int32_t fftsize, int32_t framesize)
		: SampleRate_(samplerate), FftSize_(fftsize), FrameSize_(framesize), FftBins_(fftsize / 2 + 1){}

	~VPANS();

	int32_t Initialize(int32_t* param);

	void Process(InterParamters* Param, VP_CPX* OutData, const float* Win, const int32_t WinSize);

	void Release();

private:

	void NoiseStateCalculate(InterParamters* Params);
	bool NoiseStateFlag;
	bool NoiseStateFlagIm;
	int32_t NoiseStateCounter;
	float* PnoiseBins;

	void ThreshStateCalculate(InterParamters* Params);
	bool ThreshStateFlag;
	int32_t ThreshStateCounter;
	float ThreshOut;

	void NoiseUpdate(InterParamters* Params);
	float* S;
	float* Smin;
	float* Stmp;
	int32_t nb_adapt;
	int32_t min_count;
	int32_t min_count_low;
	float* PnoiseEstimate;

	void NoiseSuppression(InterParamters* Params);
	bool zero_flag;

	void ZeroStateCalculate(InterParamters* Params);
	bool zero_flag_stable;
	int32_t ZeroCounter;

	void PrioriSNRCalcu();

	void ProbFrameCalcu(InterParamters* Params);

	void FinalGainCalcu(InterParamters* Params);

	void EchoFinalCheck(InterParamters* Params);

	void NoiseFinalCheck(InterParamters* Params);

	void FinalCheck(InterParamters* Params);

	void GetOutPut(InterParamters* Params, VP_CPX* OutData);

	int32_t IndexCheckProtect;
	int32_t IndexCheckEcho;
	int32_t IndexCheckEchoIm;

	int32_t IndexLowFreq;
	int32_t IndexHighFreq;

	int32_t NLPBins;
	float* Pnoise_nlp;
	float* Psig_nlp;
	float* WinSmooth_nlp;
	int32_t NumWinSmooth;
	float* Gain_nlp;

	float* PriorProb;
	float* PsigOld;
	float* PriorProbSmth;

	int32_t IndexLowBins;
	int32_t IndexMidBins;

	int32_t PrefStart_n;
	int32_t PrefEnd_n;
	int32_t PrefStart_n2;
	int32_t PrefEnd_n2;
	int32_t PrefStart_e;
	int32_t PrefEnd_e;
	int32_t ProbCounter;

	int32_t FinalIbeg;
	int32_t FinalIstart;
	int32_t FinalIend;

	static const float TableOUT[];
	static const float NOISETABLE[];
	static const float NOISEFLOOR[];

	uint16_t IndexNoise;
	uint16_t IndexNoise2;

	void FilterComputer(float* ps, float* mel, int32_t len);
	void FilterComputerInv(InterParamters* Params, float* mel, float* ps, int32_t len);
};