#include "vp_ns.h"
#include "vp_config.h"
#include <math.h>

VPANS::~VPANS()
{
}

int32_t VPANS::Initialize(int32_t* param)
{
	NLPBins = Min(FftBins_, (int32_t)(10000.0f * (float)FftSize_ / (float)SampleRate_ + 0.5f);

	IndexNoise = 13471;
	IndexNoise2 = 31573;

	nb_adapt = 0;
	min_count = 0;
	min_count_low = 0;

	zero_flag_stable = false;
	ZeroCounter = 0;
	zero_flag = false;
	NoiseStateFlag = false;
	NoiseStateFlagIm = false;
	NoiseStateCounter = 0;

	ThreshStateFlag = false;
	ThreshStateCounter = 0;
	ThreshOut = GAIN_THRESH_UP;

	S = new float[NLPBins];
	VectorClear(Smin, NLPBins);

	Stmp = new float[NLPBins];
	VectorClear(Stmp, NLPBins);

	Gain_nlp = new float[NLPBins];
	VectorClear(Gain_nlp, NLPBins);

	PnoiseBins = new float[NLPBins];
	VectorClear(PnoiseBins, NLPBins);

	Pnoise_nlp = new float[NLPBins];
	VectorClear(Pnoise_nlp, NLPBins);

	Psig_nlp = new float[NLPBins];
	VectorClear(Psig_nlp, NLPBins);

	PriorProb = new float[NLPBins];
	VectorClear(PriorProb, NLPBins);

	PsigOld = new float[NLPBins];
	VectorClear(PsigOld, NLPBins);

	PriorProbSmth = new float[NLPBins];
	VectorClear(PriorProbSmth, NLPBins);

	PrefStart_e = (int32_t)(100.0f *(float)FftSize_ / (float)SampleRate_ + 0.5f);
	PrefEnd_e = (int32_t)(1000.f * (float)FftSize_ / (float)SampleRate_ + 0.5f);

	PrefStart_n = (int32_t)(1800.0f * (float)FftSize_ / (float)SampleRate_ + 0.5f);
	PrefEnd_n = (int32_t)(4000.0f * (float)FftSize_ / (float)SampleRate_ + 0.5f);

	PrefStart_n2 = (int32_t)(4000.0f * (float)FftSize_ / (float)SampleRate_ + 0.5f);
	PrefEnd_n2 = (int32_t)(6500.0f * (float)FftSize_ / (float)SampleRate_ + 0.5f);

	IndexLowFreq = (int32_t)(700.0f * (float)FftSize_ / (float)SampleRate_ + 0.5f);
	IndexHighFreq = (int32_t)(4500.0f * (float)FftSize_ / (float)SampleRate_ + 0.5f);

	IndexCheckEcho = (int32_t)(700.0f * (float)FftSize_ / (float)SampleRate_ + 0.5f);
	IndexCheckEchoIm = (int32_t)(4500.0f * (float)FftSize_ / (float)SampleRate_ + 0.5f);

	IndexCheckProtect = (int32_t)(400.0f * (float)FftSize_ / (float)SampleRate_ + 0.5f);

	IndexLowBins = (int32_t)(650.f * (float)FftSize_ / (float)SampleRate_ + 0.5f);

	IndexMidBins = (int32_t)(3000.0f * (float)FftSize_ / (float)SampleRate_ + 0.5f);

	ProbCounter = 0;

	NumWinSmooth = Max(1, (int32_t)(100.0f * (float)FftSize_ / (float)SampleRate_ + 0.5f));
	WinSmooth_nlp = new float[2 * NumWinSmooth + 1];
	for (int m = 0; m < 2 * NumWinSmooth + 1; ++m) {
		WinSmooth_nlp[m] = sqrtf(0.5f * (1.0f - cosf(2.0 * M_PI * (m + 1) / (2.0f * NumWinSmooth + 2.0f))));
	}

	PnoiseEstimate = new float[NLPBins];
	VectorClear(PnoiseEstimate, NLPBins);

	FinalIbeg = (int32_t)(500.0f * (float)FftSize_ / (float)SampleRate_ + 0.5f);
	FinalIstart = (int32_t)(1500.0f * (float)FftSize_ / (float)SampleRate_ + 0.5f);
	FinalIend = Min(FftBins_, (int32_t)(4000.0f * (float)FftSize_ / (float)SampleRate_ + 0.5f));

	return 1;
}

void VPANS::Release()
{
	delete[] PnoiseEstimate;
	delete[] PnoiseBins;
	delete[] Pnoise_nlp;
	delete[] Psig_nlp;
	delete[] PsigOld;
	delete[]  PriorProb;
	delete[] PriorProbSmth;
	delete[] Gain_nlp;
	delete[] WinSmooth_nlp;

	delete[] S;
	delete[] Smin;
	delete[] Stmp;
}

void VPANS::Process(InterParamters* Param, VP_CPX* OutData, const float* Win, const int32_t WinSize)
{
	

	ThreshStateCalculate(Param);

	NoiseUpdate(Param);

	NoiseStateCalculate(Param);

	NoiseSuppression(Param);

	GetOutPut(Param, OutData);

	ZeroStateCalculate(Param);

	Param->ZeroFlag = zero_flag;
	Param->Pnoise = PnoiseEstimate;
}

void VPANS::NoiseUpdate(InterParamters* Param)
{
	float* Psig = Param->Psig;

	min_count++;
	min_count_low++;

	if (nb_adapt < 20000) {
		nb_adapt++;
	}

	float FftSize_f = 0.1f * (float)FftSize_;
	float ss = 0.0f;
	float ss2 = 0.0f;
	float ss3 = 0.0f;
	for (int32_t i = 0; i < NLPBins; ++i) {
		if (Psig[i] > NOISEFLOOR[i] * FftSize_f) {
			if (Psig[i] < 3.0f * S[i]) {
				S[i] = ss * S[i] + (1.0f - ss) * Psig[i];
			}
			else if(Psig[i] < 10.0f * S[i]) {
				S[i] = ss2 * S[i] + (1.0f - ss2) * Psig[i];
			}
			else {
				S[i] = ss3 * S[i] + (1.0f - ss3) * Psig[i];					
			}
		} 
	}

	int32_t min_range;
	if (nb_adapt < 100) {
		min_range = 15;
	}
	else if(nb_adapt < 500) {
		min_range = 100;
	}
	else {
		min_range = 200;
	}

	float Gate_d0 = 1.2f;
	float Gate_d1 = 1.5f;
	float Gate_u = 2.5f;

	/* low band */
	{
		if (min_count_low < min_range) {
			for (int32_t i = 0; i < IndexLowBins; ++i) {
				Smin[i] = Min(Smin[i], S[i]);
				Stmp[i] = Min(Stmp[i], S[i]);
			}
		}
		else {
			min_count_low = 0;
			for (int32_t i = 0; i < IndexLowBins; ++i) {
				Smin[i] = Min(Stmp[i], S[i]);
				Stmp[i] = S[i];
			}
		}

		if (!NoiseStateFlagIm) {
			for (int32_t i = 0; i < IndexLowBins; ++i) {
				if (S[i] < Gate_u * Smin[i] || Psig[i] < Gate_u * PnoiseEstimate[i]) {
					if (Psig[i] > PnoiseEstimate[i]) {
						PnoiseEstimate[i] = PnoiseEstimate[i] * 0.85f + Psig[i] * 0.15f;
					}
					else {
						PnoiseEstimate[i] = PnoiseEstimate[i] * 0.95f  + Psig[i] * 0.05f;
					}
				}
			}
		}
		else {
			for (int32_t i = 0; i < IndexLowBins; ++i) {
				if (S[i] < Gate_d0 * Smin[i] || Psig[i] < Gate_d0 * PnoiseEstimate[i]) {
					if (Psig[i] > PnoiseEstimate[i]) {
						PnoiseEstimate[i] = PnoiseEstimate[i] * 0.85f + Psig[i] * 0.15f;
					}
					else {
						PnoiseEstimate[i] = PnoiseEstimate[i] * 0.95f  + Psig[i] * 0.05f;
					}
				}
			}
		}
	}


	{
		if (min_count_low < min_range) {
			for (int32_t i = IndexLowBins; i < NLPBins; ++i) {
				Smin[i] = Min(Smin[i], S[i]);
				Stmp[i] = Min(Stmp[i], S[i]);
			}
		}
		else {
			min_count_low = 0;
			for (int32_t i = IndexLowBins; i < NLPBins; ++i) {
				Smin[i] = Min(Stmp[i], S[i]);
				Stmp[i] = S[i];
			}
		}

		if (!NoiseStateFlagIm) {
			for (int32_t i = IndexLowBins; i < NLPBins; ++i) {
				if (S[i] < Gate_u * Smin[i] || Psig[i] < Gate_u * PnoiseEstimate[i]) {
					if (Psig[i] > PnoiseEstimate[i]) {
						PnoiseEstimate[i] = PnoiseEstimate[i] * 0.85f + Psig[i] * 0.15f;
					}
					else {
						PnoiseEstimate[i] = PnoiseEstimate[i] * 0.95f  + Psig[i] * 0.05f;
					}
				}
			}
		}
		else {
			for (int32_t i = IndexLowBins; i < NLPBins; ++i) {
				if (S[i] < Gate_d1 * Smin[i] || Psig[i] < Gate_d1 * PnoiseEstimate[i]) {
					if (Psig[i] > PnoiseEstimate[i]) {
						PnoiseEstimate[i] = PnoiseEstimate[i] * 0.85f + Psig[i] * 0.15f;
					}
					else {
						PnoiseEstimate[i] = PnoiseEstimate[i] * 0.95f  + Psig[i] * 0.05f;
					}
				}
			}
		}
	}

	if (nb_adapt < 5) {
		for (int32_t i = 0; i < NLPBins; ++i) {
			Smin[i] = Max(Smin[i], Psig[i]);
			Stmp[i] = Smin[i];

			PnoiseEstimate[i] = Smin[i];
		}
	}

	FftSize_f = 0.01f * (float)FftSize_;
	int32_t IndexFreq1 = (int32_t)(1000.0f * (float)FftSize_ / (float)SampleRate_ + 0.5f);
	int32_t IndexFreq2 = (int32_t)(4500.0f * (float)FftSize_ / (float)SampleRate_ + 0.5f);
	if (zero_flag_stable) {
		for (int32_t i = 0; i < IndexFreq1; ++i) {
			PnoiseEstimate[i] = Max(PnoiseEstimate[i], NOISEFLOOR[i] * FftSize_f);
		}
		for (int32_t i = IndexFreq1; i < IndexFreq2; ++i) {
			PnoiseEstimate[i] = Max(PnoiseEstimate[i], NOISEFLOOR[i] * FftSize_f * 0.1f);
		}
	}
	else {
		for (int32_t i = 0; i < IndexFreq1; ++i) {
			PnoiseEstimate[i] = Max(PnoiseEstimate[i], NOISEFLOOR[i] * FftSize_f * 5.0f);
		}
		for (int32_t i = IndexFreq1; i < IndexFreq2; ++i) {
			PnoiseEstimate[i] = Max(PnoiseEstimate[i], NOISEFLOOR[i] * FftSize_f * 0.5f);
		}
	}
}

void VPANS::NoiseSuppression(InterParamters* Param) 
{
	for (int32_t m = 0; m < NLPBins; ++m) {
		PnoiseBins[m] = Min(Param->Precho[m] + PnoiseEstimate[m], 1.5f * Param->PsigOut[m]);
	}

	/* Transform */
	FilterComputer(PnoiseBins, Pnoise_nlp, NLPBins);
	FilterComputer(Param->PsigOut, Psig_nlp, NLPBins);

	/* Compute a priori SNR */
	PrioriSNRCalcu();

	/* Calculate the speech probability of whole frame */
	ProbFrameCalcu(Param);

	/* Computer gain */
	FinalGainCalcu(Param);

	EchoFinalCheck(Param);

	NoiseFinalCheck(Param);

	FinalCheck(Param);

	FilterComputerInv(Param, Gain_nlp, Param->GainNS, NLPBins);
}

void VPANS::PrioriSNRCalcu()
{
	for (int32_t i = 1; i < NLPBins; ++i) {
		/* Total noise estimate including residual echo and background noise */
		float Pnoise_i = Pnoise_nlp[i] + 1.0e-10f;

		/* Computering the updating paramter */
		float tmp = PsigOld[i] / (PsigOld[i] + Pnoise_i);
		float gamma = 0.19f + 0.8f * tmp * tmp;

		/* A posteriori SNR */
		float Post = Max(0.01f, Min(Psig_nlp[i] / Pnoise_i - 1.0f, 100.0f));

		/* DD: A priori SNR update */
		PriorProb[i] = gamma * Post + (1.0 - gamma) * Min(PsigOld[i] / Pnoise_i, 100.0f);

		/* A prior SNR smooth */
		PriorProbSmth[i] = 0.7f * PriorProbSmth[i] + 0.3f * PriorProb[i];
	}
}

void VPANS::ProbFrameCalcu(InterParamters* Params)
{
	float ProbFrame;

	float PriorMean = 0;
	for (int32_t i = PrefStart_e; i < PrefEnd_e; ++i) {
		PriorMean += PriorProbSmth[i];
	}
	PriorMean /= (PrefEnd_e - PrefStart_e);
	ProbFrame = 0.1f + 0.9f * PriorMean / (PriorMean + 0.5f);

	if (Params->EchoStateFlagIm) {
		PriorMean = 0;
		for (int i = PrefStart_n; i < PrefEnd_n; ++i) {
			PriorMean += PriorProbSmth[i];
		}
		PriorMean /= (PrefEnd_n - PrefStart_n);
		float ProbFrame1 = 0.1f + 0.9f * PriorMean / (PriorMean + 0.5f);

		PriorMean = 0;
		for (int i = PrefStart_n2; i < PrefEnd_n2; ++i) {
			PriorMean += PriorProbSmth[i];
		}
		PriorMean /= (PrefEnd_n2 - PrefStart_n2);
		float ProbFrame2 = 0.1f + 0.9f * PriorMean / (PriorMean + 0.5f);

		if (NoiseStateFlag) {
			float ProbSelect = Max(ProbFrame1, ProbFrame2);
			float ProbMax = Max(ProbSelect, ProbFrame);
			float ProbMin = Min(ProbSelect, ProbFrame);

			if (Params->ProbFrame < 0.5f) {
				ProbFrame = Params->ProbFrame * ProbMin + (1.0f - Params->ProbFrame) * ProbMax;
			}
			else {
				ProbFrame = Params->ProbFrame * ProbMax + (1.0f - Params->ProbFrame) * ProbMin;
			}
			ProbCounter = 0;
		}
		else {
			float ProbSelect = Max(ProbFrame1, ProbFrame2);
			float ProbMax = Max(ProbSelect, ProbFrame);
			float ProbMin = Min(ProbSelect, ProbFrame);

			if (Params->ProbFrame > 0.7f) {
				ProbCounter++;
			}
			else {
				ProbCounter = 0;
			}

			if (ProbCounter <= 3) {
				if (Params->ProbFrame > 0.5f) {
					ProbFrame = Params->ProbFrame * ProbMin + (1.0f - Params->ProbFrame) * ProbMax;
				}
				else {
					ProbFrame = Params->ProbFrame * ProbMax + (1.0f - Params->ProbFrame) * ProbMin;
				}
			}
		}
	}
	else if (NoiseStateFlag) {
		PriorMean = 0;
		for (int i = PrefStart_n; i < PrefEnd_n; ++i) {
			PriorMean += PriorProbSmth[i];
		}
		PriorMean /= (PrefEnd_n - PrefStart_n);
		float ProbFrame1 = 0.1f + 0.9f * PriorMean / (PriorMean + 0.5f);
		ProbFrame = Max(ProbFrame1, ProbFrame);
	}

	Params->ProbFrame = ProbFrame;
}


// 划分不同频带做降噪
void VPANS::FilterComputer(float* ps, float* mel, int32_t len)
{
	mel[1] = ps[1];
	mel[len - 1] = ps[len - 1];

	for (int32_t i = 2; i < NumWinSmooth; ++i) {
		mel[i] = 0;
		for (int32_t m = 1 - i; m <= i - 1; ++m) {
			mel[i] += ps[i + m] * WinSmooth_nlp[m + NumWinSmooth];
		}
	}

	for (int32_t i = len - NumWinSmooth; i < len - 1; ++i) {
		mel[i] = 0;
		for (int32_t m = -(len - i - 1); m <= len - i - 1; ++m) {
			mel[i] += ps[i + m] * WinSmooth_nlp[m + NumWinSmooth];
		}
	}

	for (int32_t i = NumWinSmooth; i < len - NumWinSmooth; ++i) {
		mel[i] = 0;
		for (int32_t m = -NumWinSmooth; m <= NumWinSmooth; ++m) {
			mel[i] += ps[i + m] * WinSmooth_nlp[m + NumWinSmooth];
		}
	}
}

void VPANS::FilterComputerInv(InterParamters* Params, float* mel, float* ps, int32_t len)
{
	if (!zero_flag_stable) {
		for (int32_t i = 0; i < FftBins_; ++i) {
			ps[i] = 0.0f;
		}
	}
	else {
		if (len < FftBins_) {
			if (!Params->EchoStateFlagIm) {
				for (int32_t i = 0; i < len; ++i) {
					ps[i] = mel[i];
				}

				for (int32_t i = len; i < FftBins_; ++i) {
					ps[i] = 0.0f;
				}
			}
			else {
				float mean = 0.0f;
				for (int32_t i = 0; i < len; ++i) {
					ps[i] = mel[i];
					mean += mel[i];
				}
				mean /= len;

				int32_t Istart = Max(0, len - (int32_t)(800.0f * (float)FftSize_ / (float)SampleRate_ + 0.5f));
				float mean2 = 0.0f;
				for (int32_t i = Istart; i < len; ++i) {
					mean2 += mel[i];
				}
				mean2 /= (len - Istart);

				mean = (mean < mean2) ? mean : mean2;

				for (int32_t i = len; i < FftBins_; ++i) {
					ps[i] = mean;
				}

			}
		}
		else {
			for (int32_t i = 0; i < len; ++i) {
				ps[i] = mel[i];
			}
		}
	}
}