#include <math.h>
#include "vp_aec.h"
#include "vp_math.h"
#include "vp_config.h"

VPAEC::~VPAEC()
{
}

int32_t VPAEC::Initialize(int32_t* param)
{
	PResidualEcho = new float[FftBins_];
	VectorClear(PResidualEcho, FftBins_);

	EC_REF_ALIGN_MAX = param[6];
	EC_REF_ALIGN_MAX = Max(15, EC_REF_ALIGN_MAX);
	EC_REF_ALIGN_MAX = Min(100, EC_REF_ALIGN_MAX);

	ThreshBand = 0.5f;
	ThreshBins = 0.7f;

	/* alignment */
	{
		MICNUM = 4;
		SEARCHNUM = EC_REF_ALIGN_MAX;
		REFNUM = SEARCHNUM + MICNUM - 1;
		SMOOTHBINS = (int32_t)(100.0f * (float)FftSize_ / (float)SampleRate_ + 0.5f);

		BinStart = (int32_t)(300.0f * (float)FftSize_ / (float)SampleRate_ + 0.5f);
		BinStart = (BinStart > (SMOOTHBINS + 1)) ? BinStart : (SMOOTHBINS + 1);

		BinStep = (int32_t)(150.0f * (float)FftSize_ / (float)SampleRate_ + 0.5f);
		int32_t BinTmp = (FftBins_ - BinStart - SMOOTHBINS - 1) / 31;
		BinStep = (BinStep < BinTmp) ? BinStep : BinTmp;

		WinSmooth = new float[2 * SMOOTHBINS + 1];
		float WinSum = 0.0f;
		for (int32_t m = 0; m < 2 * SMOOTHBINS + 1; ++m) {
			WinSmooth[m] = sqrtf(0.5f * (1.0f - cosf(2.0 * M_PI * (m + 1) / (2 * SMOOTHBINS + 2))));
			WinSum += WinSmooth[m];
 		}

		for (int32_t m = 0; m < 2 * SMOOTHBINS + 1; ++m) {
			WinSmooth[m] /= WinSum;
		}

		IndexMic = 0;
		IndexRef = 0;

		MicDataBuffer = new VP_CPX[MICNUM * FftBins_];
		VectorClear(MicDataBuffer, MICNUM * FftBins_);

		MicDataBufferAbs2 = new float[MICNUM * FftBins_];
		VectorClear(MicDataBufferAbs2, MICNUM * FftBins_);

		RefDataBuffer = new VP_CPX[REFNUM * FftBins_];
		VectorClear(RefDataBuffer, REFNUM * FftBins_);

		RefDataBufferAbs2 = new float[REFNUM * FftBins_];
		VectorClear(RefDataBufferAbs2, REFNUM * FftBins_);

		MicDataAbs2 = new float[FftBins_];
		VectorClear(MicDataAbs2, FftBins_);

		RefDataAbs2 = new float[FftBins_];
		VectorClear(RefDataAbs2, FftBins_);

		MicFeature = new Feature[MICNUM];
		for (int m = 0; m < MICNUM; ++m) {
			MicFeature[m].valid = false;
			MicFeature[m].fe = 0;
		}

		RefFeature = new Feature[REFNUM];
		for (int m = 0; m < SEARCHNUM; ++m) {
			RefFeature[m].valid = false;
			RefFeature[m].fe = 0;
		}

		FeatureCurrent = new float[SEARCHNUM + 2];
		VectorClear(FeatureCurrent, SEARCHNUM + 2);

		FeatureSmooth = new float[SEARCHNUM + 2];
		for (int m = 0; m < SEARCHNUM; ++m) {
			FeatureSmooth[m] = 0.7f;
		}

		FeatureCopy = new float[SEARCHNUM + 2];
		VectorClear(FeatureCopy, SEARCHNUM + 2);

		SmoothCounter = 0;
		IndexAlignment = 0;
	}

	/* linear ec */
	{
		NumStage_a = 3;
		DetectionFlag = new float[FftBins_];
		VectorClear(DetectionFlag, FftBins_);

		Wa = new VP_CPX[SEARCHNUM * NumStage_a * FftBins_];
		VectorClear(Wa, SEARCHNUM * NumStage_a * FftBins_);

		Pmr = new VP_CPX[SEARCHNUM * NumStage_a * FftBins_];
		VectorClear(Pmr, SEARCHNUM * NumStage_a * FftBins_);

		Prr = new float[SEARCHNUM * NumStage_a * FftBins_];
		VectorClear(Prr, SEARCHNUM * NumStage_a * FftBins_);

		IndexMic4EC = new int[3];
		IndexRef4EC = new int[NumStage_a];

		CorrTmp1 = new float[FftBins_];
		VectorClear(CorrTmp1, FftBins_);

		CorrTmp2 = new float[FftBins_];
		VectorClear(CorrTmp2, FftBins_);

		CorrTmpc = new VP_CPX[FftBins_];
		VectorClear(CorrTmpc, FftBins_);

		DetectionBins_s = (int32_t)(500.0f *(float)FftSize_ / (float)SampleRate_ + 0.5f);
		DetectionBins_e = (int32_t)(2000.0f * (float)FftSize_ / (float)SampleRate_ + 0.5f);

		EchoEst = new VP_CPX[FftBins_];
		VectorClear(EchoEst, FftBins_);

		EchoEstAbs2 = new float[FftBins_];
		VectorClear(EchoEstAbs2, FftBins_);

		EchoEst_bk = new VP_CPX[FftBins_];
		VectorClear(EchoEst_bk, FftBins_);

		OutData_bk = new VP_CPX[FftBins_];
		VectorClear(OutData_bk, FftBins_);

		Pecho4ES = 0.0f;
		Psig4ES = 0.0f;
		Counter4ES = 0;

		Index4ES_s = (int)(600.0f *(float)FftSize_ / SampleRate_ + 0.5f);
		Index4ES_e = (int)(3900.0f * (float)FftSize_ / SampleRate_ + 0.5f);

		IndexLowFreq = (int)(600.0f *(float)FftSize_ / SampleRate_ + 0.5f);
		IndexMidFreq = Min((int)(3000.0f * (float)FftSize_ / (float)SampleRate_ + 0.5f), FftBins_ - 5);
	}

	/* search ec */
	{
		SearchFrameValid = new int32_t[SEARCHNUM];
		IndexAlignment = 0;

		SearchStart = (int32_t)(800.0f * (float)FftSize_ / (float)SampleRate_ + 0.5f);
		SearchEnd = (int32_t)(6000.0f * (float)FftSize_ / (float)SampleRate_ + 0.5f);

		SearchNum = SearchEnd - SearchStart + 1;
		SearchNum = (SearchNum + 1) / 2;

		NumStage_s = 2;

		Wa_s = new VP_CPX[SEARCHNUM * NumStage_s * SearchNum];
		VectorClear(Wa_s, SEARCHNUM * NumStage_s * SearchNum);

		Pmr_s = new VP_CPX[SEARCHNUM * NumStage_s * SearchNum];
		VectorClear(Pmr_s, SEARCHNUM * NumStage_s * SearchNum);

		Prr_s = new float[SEARCHNUM * NumStage_s * SearchNum];
		VectorClear(Prr_s, SEARCHNUM * NumStage_s * SearchNum);

		CancelOutSumCurr = new float[SEARCHNUM];
		VectorClear(CancelOutSumCurr, SEARCHNUM);

		CancelOutSumSmth = new float[SEARCHNUM];
		VectorClear(CancelOutSumSmth, SEARCHNUM);
	}

	/* NLP */
	{
		PostBins_s = (int32_t)(600.0f * (float)FftSize_ / (float)SampleRate_ + 0.5f);
		PostBins_e = (int32_t)(3000.0f * (float)FftSize_ / (float)SampleRate_ + 0.5f);
		PostState = 0;
		PostStateCounter = 0;

		PechoSmooth = 1.0f;
		PsigSmooth = 1.0f;
		SEROut = 1.0f;

		Ps = new float[FftBins_];
		VectorClear(Ps, FftBins_);

		Pe = new float[FftBins_];
		VectorClear(Pe, FftBins_);
		
		Ps0 = new float[FftBins_];
		VectorClear(Ps0, FftBins_);

		Pe0 = new float[FftBins_];
		VectorClear(Pe0, FftBins_);

		Gain0 = new float[FftBins_];
		VectorClear(Gain0, FftBins_);

		EchoStateCounter = 0;

		PResidualEcho0 = new float[FftBins_];
		VectorClear(PResidualEcho0, FftBins_);

	}

	return 1;
}

void VPAEC::Release()
{
	delete[] PResidualEcho0;

	if (AECValid != 1) {
		return;
	}

	delete[] WinSmooth;
	delete[] MicDataBuffer;
	delete[] RefDataBuffer;
	delete[] MicDataBufferAbs2;
	delete[] RefDataBufferAbs2;

	delete[] MicDataAbs2;
	delete[] RefDataAbs2;

	delete[] MicFeature;
	delete[] RefFeature;

	delete[] FeatureCurrent;
	delete[] FeatureSmooth;
	delete[] FeatureCopy;

	delete[] IndexMic4EC;
	delete[] IndexRef4EC;

	delete[] CorrTmp1;
	delete[] CorrTmp2;
	delete[] CorrTmpc;

	delete[] EchoEst;
	delete[] EchoEst_bk;
	delete[] EchoEstAbs2;

	delete[] Wa_s;
	delete[] Pmr_s;
	delete[] Prr_s;
	delete[] SearchFrameValid;
	delete[] CancelOutSumCurr;
	delete[] CancelOutSumSmth;

	delete[] Ps;
	delete[] Pe;
	delete[] Ps0;
	delete[] Pe0;
	delete[] Gain0;

	delete[] PResidualEcho0;
}

void VPAEC::Process(InterParamters* Params, const VP_CPX * MicData, const VP_CPX* RefData, VP_CPX* OutData)
{
	if (AECValid != 1) {
		Params->EchoStateFlag = true;
		Params->EchoStateFlagIm = true;

		for (int m =0; m < FftBins_; ++m) {
			OutData[m] = VP_CPX(MicData[m].re(), MicData[m].im());
			Params->GainEC[m] = 1.0f;

			Params->PsigOut[m] = OutData[m].re() * OutData[m].re() + OutData[m].im() * OutData[m].im();

			Params->PsigIn[m] = Params->PsigOut[m];
		}

		Params->Psig[1] = 0.5f * Params->PsigOut[1] + 0.3f * Params->PsigOut[2] + 0.2f * Params->PsigOut[3];
		Params->Psig[2] = 0.2f * Params->PsigOut[1] + 0.3f * Params->PsigOut[2] + 0.3f * Params->PsigOut[3] + 0.2f * Params->PsigOut[4];
		Params->Psig[FftBins_ - 1] = 0.5 * Params->PsigOut[FftBins_ - 1] + 0.3f * Params->PsigOut[FftBins_ - 2] + 0.2f * Params->PsigOut[FftBins_ - 3];
		Params->Psig[FftBins_ - 2] = 0.2f * Params->PsigOut[FftBins_ - 1] + 0.3f * Params->PsigOut[FftBins_ -2] + 0.3f * Params->PsigOut[FftBins_ - 3] + 0.2f * Params->PsigOut[FftBins_ - 4];
		for (int m = 3; m < FftBins_ - 2; ++m) {
			Params->Psig[m] = Params->PsigOut[m - 2] * 0.15f + Params->PsigOut[m - 1] * 0.2f
				+ Params->PsigOut[m] * 0.3f + Params->PsigOut[m + 1] * 0.2f + Params->PsigOut[m + 2] * 0.15f;
		}

		Params->Precho = PResidualEcho;

		return;
	}

	IndexMic = GetMicIndex(-1);
	IndexRef = GetRefIndex(-1);

	for (int m = 1; m < FftBins_; ++m) {
		MicDataAbs2[m] = MicData[m].re() * MicData[m].re() + MicData[m].im() * MicData[m].im();
		RefDataAbs2[m] = RefData[m].re() * RefData[m].re() + RefData[m].im() * RefData[m].im();

		MicDataBuffer[m * MICNUM + IndexMic] = MicData[m];
		RefDataBuffer[m * REFNUM + IndexRef] = RefData[m];

		MicDataBufferAbs2[m * MICNUM + IndexMic] = MicDataAbs2[m];
		RefDataBufferAbs2[m * REFNUM + IndexRef] = RefDataAbs2[m];

		Params->PsigIn[m] = MicDataAbs2[m];
	}

	for (int m = 0; m < 3; ++m) {
		IndexMic4EC[m] = GetMicIndex(m);
	}

	float Bias = 0.1f * (float)FftSize_;
	MicFeature[IndexMic] = FeatureCaculate(MicDataAbs2, -Bias);
	RefFeature[IndexRef] = FeatureCaculate(RefDataAbs2, Bias);
 
	FeatureCurrCalcu();

	FeatureSmthCalcu();

	EchoCanceling(Params, OutData);

	EchoStateCalculate(Params);

	EchoSuppression(Params);
}

void VPAEC::EchoStateCalculate(InterParamters* Params)
{
	float Pecho = 0.0f;
	float Psig = 0.0f;

	for (int m = Index4ES_s; m <= Index4ES_e; ++m) {
		Pecho += EchoEstAbs2[m];
		Psig += MicDataAbs2[m];
	}

	Pecho4ES = Pecho4ES * 0.5f + Pecho;
	Psig4ES = Psig4ES * 0.5f + Psig;

	switch (Params->EchoStateFlag)
	{
	case true:
		if (Pecho4ES > 0.1f *(Psig4ES + 1.0e-4f)) {
			Counter4ES++;
			Params->EchoStateFlagIm = false;
		}
		else {
			Counter4ES = 0;
			Params->EchoStateFlagIm = true;
		}
		if (Counter4ES >= 3) {
			Counter4ES = 0;
			Params->EchoStateFlag = false;
		}
		EchoStateCounter = 0;
		break;

	case false:
		if (Pecho4ES < 0.1f * (Psig4ES + 1.0e-4f)) {
			Counter4ES++;
			Params->EchoStateFlagIm = true;
		}
		else {
			Counter4ES = 0;
			Params->EchoStateFlagIm = false;
		}

		if(Counter4ES >= 20) {
			Counter4ES = 0;
			Params->EchoStateFlag = true;
		}

		if (Params->ProbFrame > 0.4f) {
			EchoStateCounter++;
		}
		else {
			EchoStateCounter = 0;
		}

	default:

		Params->EchoStateFlag = false;
		Params->EchoStateFlagIm = false;

		break;
	}
}

void VPAEC::FeatureCurrCalcu()
{
	FeatureCurrent[SEARCHNUM] = 0.0f;

	float fetotal = (float)(32.0f * MICNUM);

	int IndexMax = 0;
	for (int index = 0; index < SEARCHNUM; ++index) {
		int value = 0;
		for (int i = 0; i < MICNUM; ++i) {
			int ei = GetMicIndex(i);
			int ri = GetRefIndex(index + i);

			if (MicFeature[ei].valid && RefFeature[ri].valid) {
				value += 32 - Bcounter(MicFeature[ei].fe ^ RefFeature[ri].fe);
			}
		}

		FeatureCurrent[index] = ((float)value) / fetotal;

		if (FeatureCurrent[SEARCHNUM] < FeatureCurrent[index]) {
			FeatureCurrent[SEARCHNUM] = FeatureCurrent[index];
			IndexMax = index;
		}
	}

	FeatureCurrent[SEARCHNUM + 1] = (float)IndexMax;
}

void VPAEC::FeatureSmthCalcu()
{
	if (FeatureCurrent[SEARCHNUM] > 0.7f) {
		float s = 0.993f;
		if (IndexAlignment != IndexAlignSearch) {
			s = 0.98f;
		}
		if (SmoothCounter < 10) {
			SmoothCounter++;
			s = 0.9f;
		}

		FeatureSmooth[SEARCHNUM] = 0.0f;
		int IndexMax;
		for (int index = 0; index < SEARCHNUM; ++index) {
			FeatureSmooth[index] = FeatureSmooth[index] * s + FeatureCurrent[index] * (1.0f - s);

			if (FeatureSmooth[SEARCHNUM] < FeatureSmooth[index]) {
				FeatureSmooth[SEARCHNUM] = FeatureSmooth[index];
				IndexMax = index;
			}
		}

		IndexAlignment = IndexMax;
		FeatureSmooth[SEARCHNUM + 1] = (float)IndexMax;
	}
}

void VPAEC::EchoSuppression(InterParamters* Params)
{
	if (!Params->EchoStateFlag) {
		float Pecho = 0.0f;
		float Psig = 0.0f;
		for (int m = Index4ES_s; m < Index4ES_e; ++m) {
			Pecho += Min(EC_REF_ALIGN_MAX * EchoEstAbs2[m], Params->PsigOut[m]);
			Psig += Params->PsigOut[m];
		}

		PechoSmooth = PechoSmooth * 0.95f + Pecho * 0.05f;
		PsigSmooth = PsigSmooth * 0.95f + Psig * 0.05f;

		SEROut = PechoSmooth / (PsigSmooth + 1.0e-10f);
	}

	float ratio = Max(5.0f * SEROut - 0.5f, EC_MIN_GATE);
	float s0 = 0.75f;   // 快平滑
	float s0_1 = 1.0f - s0;
	float s = 0.9f;    // 慢平滑
	float s_1 = 1.0f - s;

	for (int m = 0; m < FftBins_; ++m) {
		float tmp = ratio * EchoEstAbs2[m];
		tmp = Min(tmp, EC_MIN_GATE * Params->PsigOut[m]);
		PResidualEcho0[m] = Max(PResidualEcho0[m] * s0, tmp);
		PResidualEcho[m] = Max(PResidualEcho[m] * s, tmp);
	}

	if (!Params->EchoStateFlagIm && Params->EchoStateFlag) {
		for (int m  = 0; m < FftBins_; ++m) {
			float tmp = Params->PsigOut[m];

			Ps0[m] = Ps0[m] * s0 + tmp * s0_1;
			Pe0[m] = Pe0[m] * s0 + Min(tmp, PResidualEcho0[m]) * s0_1;
			// 当前信回比
			float SER2_c = Max(0.001f, (tmp + 1.0e-5f) / (Min(tmp, PResidualEcho0[m]) + 1.0e-10f) - 1.0f);
			// 平滑信回比
			float SER2 = Min(SER2_c, Max(0.001f, (Ps0[m] + 1.0e-5f) / (Pe0[m] + 1.0e-10f)) - 1.0f);
			Gain0[m] = SER2 / (SER2 + 1.0f);

			Ps[m] = Ps[m] * s + tmp * s_1;
			Pe[m] = Pe[m] * s + Min(tmp, PResidualEcho[m]) * s_1;
			float SER_c = Max(0.001f, (tmp + 1.0e-5f) / (Min(tmp, PResidualEcho[m]) + 1.0e-10f) - 1.0f);
			float SER = Min(SER_c, Max(0.001f, (Ps[m] + 1.0e-5f) / (Pe[m] + 1.0e-10f)) - 1.0f);
			Params->GainEC[m] = SER / (SER + 1.0f);
		}
	}
	else {
		for (int m  = 0; m < FftBins_; ++m) {
			float tmp = Params->PsigOut[m];

			Ps0[m] = Ps0[m] * s0 + tmp * s0_1;
			Pe0[m] = Pe0[m] * s0 + Min(tmp, PResidualEcho0[m]) * s0_1;
			// 平滑信回比
			float SER2 = Max(0.001f, (Ps0[m] + 1.0e-5f) / (Pe0[m] + 1.0e-10f) - 1.0f);
			Gain0[m] = SER2 / (SER2 + 1.0f);

			Ps[m] = Ps[m] * s + tmp * s_1;
			Pe[m] = Pe[m] * s + Min(tmp, PResidualEcho[m]) * s_1;
			float SER =Max(0.001f, (Ps[m] + 1.0e-5f) / (Pe[m] + 1.0e-10f) - 1.0f);
			Params->GainEC[m] = SER / (SER + 1.0f);
		}
	}

	if (EC_POST_LEVEL == 0) {
		for (int m = 0; m < FftBins_; ++m) {
			Params->GainEC[m] = Gain0[m];
		}
	}
	else if (EC_POST_LEVEL == 2) {
		for (int m = 0; m < FftBins_; ++m) {
			Params->GainEC[m] = Min(Params->GainEC[m], Gain0[m]);
		}
	}
	else if (EC_POST_LEVEL == 1) {
		int PostBinsNum0 = 0;
		int PostBinsNum1 = 0;
		int PostAll = PostBins_e - PostBins_s;

		switch (PostState)
		{
		case 0:
			for (int m = PostBins_s; m < PostBins_e; ++m) {
				if (Gain0[m] < 0.1f) {
					PostBinsNum0++;
				}
				if (Params->GainEC[m] < 0.1f) {
					PostBinsNum1++;
				}
			}

			if (PostBinsNum0 > 0.5f * PostAll && PostBinsNum1 > 0.75f * PostAll) {
				PostStateCounter++;
			}
			else {
				PostStateCounter = 0;
			}


			if (PostStateCounter >= 5) {
				PostStateCounter = 0;
				PostState = 1;
			}
			break;

		case 1:
			for (int m = PostBins_s; m < PostBins_e; ++m) {
				if (Gain0[m] < 0.1f) {
					PostBinsNum0++;
				}
				if (Params->GainEC[m] < 0.1f) {
					PostBinsNum1++;
				}
			}

			if (PostBinsNum0 > 0.6f * PostAll && PostBinsNum1 > 0.4f * PostAll) {
				PostStateCounter++;
			}
			else {
				PostStateCounter = 0;
			}


			if (PostStateCounter >= 2) {
				PostStateCounter = 0;
				PostState = 0;
			}
			break;

		default:
			break;
		}

		if (PostState == 0) {
			for (int m  = 0; m < FftBins_; ++m){
				Params->GainEC[m] = Gain0[m];
			}

		}
		else {
			for (int m = 0; m < FftBins_; ++m) {
				Params->GainEC[m] = Min(Params->GainEC[m], Gain0[m]);
			}
		}

	}

	Params->Precho = PResidualEcho0;
	if (!Params->EchoStateFlagIm || (!Params->EchoStateFlag && EchoStateCounter <= 3)) {
		Params->Precho = PResidualEcho;
	}
}

int VPAEC::GetRefIndex(int i)
{
	return (IndexRef + REFNUM - i) % REFNUM;
}

int VPAEC::GetMicIndex(int i)
{
	return (IndexMic + MICNUM) % MICNUM;
}

int VPAEC::Bcounter(uint32_t n)
{
	n = (n & B_01010101) + ((n >> 1) & B_01010101);
	n = (n & B_00111001) + ((n >> 2) & B_00111001);
	n = (n & B_00001111) + ((n >> 4) & B_00001111);
	return n % 255;
}

Feature VPAEC::FeatureCaculate(const float* data2, const float bias)
{
	Feature result;
	
	uint32_t xfe = 0;
	uint32_t xbs = 1;
	float bias_act = bias * 1.0e-3f;
	for (int m = 0; m < 32; ++m) {
		int idx = BinStart + m * BinStep;
		float x = data2[idx] + bias_act;

		int idx_s = idx - SMOOTHBINS;
		float x_smooth = 0.0f;
		for (int n = 0; n <= SMOOTHBINS * 2; ++n) {
			x_smooth += data2[idx_s + n] *WinSmooth[n];
		}
		xfe = (x_smooth < x) ? (xfe + xbs) : xfe;

		xbs = xbs * 2;
	}
	result.fe = xfe;

	float rate = ((float)Bcounter(xfe)) / 32.0f;
	result.valid = (rate < 0.9f && rate > 0.1f) ? true : false;

	return result;
}

void VPAEC::DetectModule4Search()
{
	float SumUp = 0.0f;
	float SumDn = 1.0e-10f;
	for(int km = 0; km < SearchNum / 3; ++km) {
		int k = km * 2 + SearchStart;
		VP_CPX* MicDataBuffer_p = MicDataBuffer + k * MICNUM;
		VP_CPX* RefDataBuffer_p = RefDataBuffer + k * REFNUM;
		float* MicDataBuffer2_p = MicDataBufferAbs2 + k * MICNUM;
		float* RefDataBuffer2_p = RefDataBufferAbs2 + k * REFNUM;

		CorrTmp1[k] = MicDataBuffer2_p[IndexMic4EC[0]];
		CorrTmp2[k] = RefDataBuffer2_p[IndexRef4EC[0]];

		CorrTmpc[k] = VP_CPX(MicDataBuffer_p[IndexMic4EC[0]].re() * RefDataBuffer_p[IndexRef4EC[0]].re()
			+ MicDataBuffer_p[IndexMic4EC[0]].im() * RefDataBuffer_p[IndexRef4EC[0]].im(),
			-MicDataBuffer_p[IndexMic4EC[0]].re() * RefDataBuffer_p[IndexRef4EC[0]].im()
			+ MicDataBuffer_p[IndexMic4EC[0]].im() * RefDataBuffer_p[IndexRef4EC[0]].re());

		for (int m = 1; m < 3; ++m) {
			CorrTmp1[k] += MicDataBuffer2_p[IndexMic4EC[m]];
			CorrTmp2[k] += RefDataBuffer2_p[IndexRef4EC[m]];
			VP_CPX tmp = VP_CPX(MicDataBuffer_p[IndexMic4EC[m]].re() * RefDataBuffer_p[IndexRef4EC[m]].re()
				+ MicDataBuffer_p[IndexMic4EC[m]].im() * RefDataBuffer_p[IndexRef4EC[m]].im(),
				-MicDataBuffer_p[IndexMic4EC[m]].re() * RefDataBuffer_p[IndexRef4EC[m]].im()
				+ MicDataBuffer_p[IndexMic4EC[m]].im() * RefDataBuffer_p[IndexRef4EC[m]].re());
			CorrTmpc[k] = VP_CPX(CorrTmpc[k].re() + tmp.re(), CorrTmpc[k].im() + tmp.im());
		}

		CorrTmp1[k] *= CorrTmp2[k];
		CorrTmp2[k] = CorrTmpc[k].re() * CorrTmpc[k].re() + CorrTmpc[k].im() * CorrTmpc[k].im();
		SumUp += CorrTmp2[k];
		SumDn += CorrTmp1[k];
	}

	float CorrBand = SumUp / SumDn;
	if (CorrBand > 0.6f * ThreshBand) {
		for( int km = SearchNum / 3; km < SearchNum; ++km) {
			int k = km * 2 + SearchStart;
			VP_CPX* MicDataBuffer_p = MicDataBuffer + k * MICNUM;
			VP_CPX* RefDataBuffer_p = RefDataBuffer + k * REFNUM;
			float* MicDataBuffer2_p = MicDataBufferAbs2 + k * MICNUM;
			float* RefDataBuffer2_p = RefDataBufferAbs2 + k * REFNUM;

			CorrTmp1[k] = MicDataBuffer2_p[IndexMic4EC[0]];
			CorrTmp2[k] = RefDataBuffer2_p[IndexRef4EC[0]];

			CorrTmpc[k] = VP_CPX(MicDataBuffer_p[IndexMic4EC[0]].re() * RefDataBuffer_p[IndexRef4EC[0]].re()
				+ MicDataBuffer_p[IndexMic4EC[0]].im() * RefDataBuffer_p[IndexRef4EC[0]].im(),
				-MicDataBuffer_p[IndexMic4EC[0]].re() * RefDataBuffer_p[IndexRef4EC[0]].im()
				+ MicDataBuffer_p[IndexMic4EC[0]].im() * RefDataBuffer_p[IndexRef4EC[0]].re());

			for (int m = 1; m < 3; ++m) {
				CorrTmp1[k] += MicDataBuffer2_p[IndexMic4EC[m]];
				CorrTmp2[k] += RefDataBuffer2_p[IndexRef4EC[m]];
				VP_CPX tmp = VP_CPX(MicDataBuffer_p[IndexMic4EC[m]].re() * RefDataBuffer_p[IndexRef4EC[m]].re()
					+ MicDataBuffer_p[IndexMic4EC[m]].im() * RefDataBuffer_p[IndexRef4EC[m]].im(),
					-MicDataBuffer_p[IndexMic4EC[m]].re() * RefDataBuffer_p[IndexRef4EC[m]].im()
					+ MicDataBuffer_p[IndexMic4EC[m]].im() * RefDataBuffer_p[IndexRef4EC[m]].re());
				CorrTmpc[k] = VP_CPX(CorrTmpc[k].re() + tmp.re(), CorrTmpc[k].im() + tmp.im());
			}

			CorrTmp1[k] *= CorrTmp2[k];
			CorrTmp2[k] = CorrTmpc[k].re() * CorrTmpc[k].re() + CorrTmpc[k].im() * CorrTmpc[k].im();
		}

		DetectionFlag[0] = CorrBand;
		for (int km = 1; km < SearchNum; ++km) {
			int k = km * 2 + SearchStart;
			DetectionFlag[k] = (MicDataBufferAbs2[k * MICNUM + IndexMic4EC[0]] > RefDataBufferAbs2[k * REFNUM + IndexRef4EC[0]] * EC_BOUND_THRESH) ?
				0.0f : (CorrTmp2[k] / (CorrTmp1[k] + 1.0e-10f));
		}
	}
	else {
		for( int km = SearchNum / 3; km < SearchNum; ++km) {
			int k = km * 2 + SearchStart;
			VP_CPX* MicDataBuffer_p = MicDataBuffer + k * MICNUM;
			VP_CPX* RefDataBuffer_p = RefDataBuffer + k * REFNUM;
			float* MicDataBuffer2_p = MicDataBufferAbs2 + k * MICNUM;
			float* RefDataBuffer2_p = RefDataBufferAbs2 + k * REFNUM;

			CorrTmp1[k] = MicDataBuffer2_p[IndexMic4EC[0]];
			CorrTmp2[k] = RefDataBuffer2_p[IndexRef4EC[0]];

			CorrTmpc[k] = VP_CPX(MicDataBuffer_p[IndexMic4EC[0]].re() * RefDataBuffer_p[IndexRef4EC[0]].re()
				+ MicDataBuffer_p[IndexMic4EC[0]].im() * RefDataBuffer_p[IndexRef4EC[0]].im(),
				-MicDataBuffer_p[IndexMic4EC[0]].re() * RefDataBuffer_p[IndexRef4EC[0]].im()
				+ MicDataBuffer_p[IndexMic4EC[0]].im() * RefDataBuffer_p[IndexRef4EC[0]].re());

			for (int m = 1; m < 3; ++m) {
				CorrTmp1[k] += MicDataBuffer2_p[IndexMic4EC[m]];
				CorrTmp2[k] += RefDataBuffer2_p[IndexRef4EC[m]];
				VP_CPX tmp = VP_CPX(MicDataBuffer_p[IndexMic4EC[m]].re() * RefDataBuffer_p[IndexRef4EC[m]].re()
					+ MicDataBuffer_p[IndexMic4EC[m]].im() * RefDataBuffer_p[IndexRef4EC[m]].im(),
					-MicDataBuffer_p[IndexMic4EC[m]].re() * RefDataBuffer_p[IndexRef4EC[m]].im()
					+ MicDataBuffer_p[IndexMic4EC[m]].im() * RefDataBuffer_p[IndexRef4EC[m]].re());
				CorrTmpc[k] = VP_CPX(CorrTmpc[k].re() + tmp.re(), CorrTmpc[k].im() + tmp.im());
			}

			CorrTmp1[k] *= CorrTmp2[k];
			CorrTmp2[k] = CorrTmpc[k].re() * CorrTmpc[k].re() + CorrTmpc[k].im() * CorrTmpc[k].im();
		}

		
		DetectionFlag[0] = CorrBand;
		for (int km = 0; km < SearchNum; ++km) {
			int k = km * 2 + SearchStart;
			DetectionFlag[k] = (MicDataBufferAbs2[k * MICNUM + IndexMic4EC[0]] > RefDataBufferAbs2[k * REFNUM + IndexRef4EC[0]] * EC_BOUND_THRESH) ?
				0.0f : 0.83f * (CorrTmp2[k] / (CorrTmp1[k] + 1.0e-10f));
		}
	}
}

void VPAEC::EchoCanceling(InterParamters* Params, VP_CPX* OutData)
{
	EchoCancelingProcess4Search();
	EchoCancelingProcess(OutData, EchoEst, IndexAlignment);

	if (IndexAlignSearch == IndexAlignment) {
		for (int m = 1; m < FftBins_; ++m) {
			Params->PsigOut[m] = OutData[m].re() * OutData[m].re() + OutData[m].im() * OutData[m].im();
			if (MicDataAbs2[m] < Params->PsigOut[m]) {
				OutData[m] = MicDataBuffer[m * MICNUM + IndexMic4EC[0]];
				Params->PsigOut[m] = MicDataAbs2[m];
			}
		}
	}
	else {
		EchoCancelingProcess(OutData_bk, EchoEst_bk, IndexAlignSearch);
		for (int m = 0; m < FftBins_; ++m) {
			Params->PsigOut[m] = OutData[m].re() * OutData[m].re() + OutData[m].im() * OutData[m].im();
			float tmp_bk = OutData_bk[m].re() * OutData_bk[m].re() + OutData_bk[m].im() * OutData_bk[m].im();
			if (Params->PsigOut[m] > MicDataAbs2[m] && tmp_bk > MicDataAbs2[m]) {
				OutData[m] = MicDataBuffer[m * MICNUM + IndexMic4EC[0]];
				Params->PsigOut[m] = MicDataAbs2[m];
			}
			else if(tmp_bk < Params->PsigOut[m]) {
				OutData[m] = OutData_bk[m];
				EchoEst[m] = EchoEst_bk[m];
				Params->PsigOut[m] = tmp_bk;
			}
		}
	}

	Params->Psig[1] = 0.5f * Params->PsigOut[1] + 0.3f * Params->PsigOut[2] 
		+ 0.2f * Params->PsigOut[3];
	Params->Psig[2] = 0.2f * Params->PsigOut[1] + 0.3f * Params->PsigOut[2] 
		+ 0.3f * Params->PsigOut[3]+ 0.2f * Params->PsigOut[4];
	Params->Psig[FftBins_ - 1] = 0.5f * Params->PsigOut[FftBins_ - 1] + 0.3f * Params->PsigOut[FftBins_ - 2] 
		+ 0.2f * Params->PsigOut[FftBins_ - 3];
	Params->Psig[FftBins_ - 2] = 0.3f * Params->PsigOut[FftBins_ - 1] + 0.3f * Params->PsigOut[FftBins_ - 2] 
		+ 0.2f * Params->PsigOut[FftBins_ - 3] + 0.3f * Params->PsigOut[FftBins_ - 4];
	
	for (int m = 3; m < FftBins_ - 2; ++m) {
		Params->Psig[m] = 0.15f * Params->PsigOut[m - 2] + 0.2f * Params->PsigOut[m - 1] 
		+ 0.3f * Params->PsigOut[m] + 0.2f * Params->PsigOut[m + 1]+ 0.15f * Params->PsigOut[m + 2];
	}

	for (int m = 1; m < FftBins_; ++m) {
		EchoEstAbs2[m] = Min(EchoEst[m].re() * EchoEst[m].re() + EchoEst[m].im() * EchoEst[m].im(), MicDataAbs2[m]); 
	}
}


void VPAEC::DetectModule()
{
	float SumUp = 0.0f;
	float SumDn = 1.0e-10f;
	for( int k = DetectionBins_s; k < DetectionBins_e; ++k) {
		VP_CPX* MicDataBuffer_p = MicDataBuffer + k * MICNUM;
		VP_CPX* RefDataBuffer_p = RefDataBuffer + k * REFNUM;
		float* MicDataBuffer2_p = MicDataBufferAbs2 + k * MICNUM;
		float* RefDataBuffer2_p = RefDataBufferAbs2 + k * REFNUM;

		CorrTmp1[k] = MicDataBuffer2_p[IndexMic4EC[0]];
		CorrTmp2[k] = RefDataBuffer2_p[IndexRef4EC[0]];

		CorrTmpc[k] = VP_CPX(MicDataBuffer_p[IndexMic4EC[0]].re() * RefDataBuffer_p[IndexRef4EC[0]].re()
			+ MicDataBuffer_p[IndexMic4EC[0]].im() * RefDataBuffer_p[IndexRef4EC[0]].im(),
			-MicDataBuffer_p[IndexMic4EC[0]].re() * RefDataBuffer_p[IndexRef4EC[0]].im()
			+ MicDataBuffer_p[IndexMic4EC[0]].im() * RefDataBuffer_p[IndexRef4EC[0]].re());

		for (int m = 1; m < 3; ++m) {
			CorrTmp1[k] += MicDataBuffer2_p[IndexMic4EC[m]];
			CorrTmp2[k] += RefDataBuffer2_p[IndexRef4EC[m]];
			VP_CPX tmp = VP_CPX(MicDataBuffer_p[IndexMic4EC[m]].re() * RefDataBuffer_p[IndexRef4EC[m]].re()
				+ MicDataBuffer_p[IndexMic4EC[m]].im() * RefDataBuffer_p[IndexRef4EC[m]].im(),
				-MicDataBuffer_p[IndexMic4EC[m]].re() * RefDataBuffer_p[IndexRef4EC[m]].im()
				+ MicDataBuffer_p[IndexMic4EC[m]].im() * RefDataBuffer_p[IndexRef4EC[m]].re());
			CorrTmpc[k] = VP_CPX(CorrTmpc[k].re() + tmp.re(), CorrTmpc[k].im() + tmp.im());
		}

		CorrTmp1[k] *= CorrTmp2[k];
		CorrTmp2[k] = CorrTmpc[k].re() * CorrTmpc[k].re() + CorrTmpc[k].im() * CorrTmpc[k].im();
		SumUp += CorrTmp2[k];
		SumDn += CorrTmp1[k];
	}

	float CorrBand = SumUp / SumDn;
	if (CorrBand > ThreshBand) {
		for( int k = 1; k < DetectionBins_s; ++k) {
			VP_CPX* MicDataBuffer_p = MicDataBuffer + k * MICNUM;
			VP_CPX* RefDataBuffer_p = RefDataBuffer + k * REFNUM;
			float* MicDataBuffer2_p = MicDataBufferAbs2 + k * MICNUM;
			float* RefDataBuffer2_p = RefDataBufferAbs2 + k * REFNUM;

			CorrTmp1[k] = MicDataBuffer2_p[IndexMic4EC[0]];
			CorrTmp2[k] = RefDataBuffer2_p[IndexRef4EC[0]];

			CorrTmpc[k] = VP_CPX(MicDataBuffer_p[IndexMic4EC[0]].re() * RefDataBuffer_p[IndexRef4EC[0]].re()
				+ MicDataBuffer_p[IndexMic4EC[0]].im() * RefDataBuffer_p[IndexRef4EC[0]].im(),
				-MicDataBuffer_p[IndexMic4EC[0]].re() * RefDataBuffer_p[IndexRef4EC[0]].im()
				+ MicDataBuffer_p[IndexMic4EC[0]].im() * RefDataBuffer_p[IndexRef4EC[0]].re());

			for (int m = 1; m < 3; ++m) {
				CorrTmp1[k] += MicDataBuffer2_p[IndexMic4EC[m]];
				CorrTmp2[k] += RefDataBuffer2_p[IndexRef4EC[m]];
				VP_CPX tmp = VP_CPX(MicDataBuffer_p[IndexMic4EC[m]].re() * RefDataBuffer_p[IndexRef4EC[m]].re()
					+ MicDataBuffer_p[IndexMic4EC[m]].im() * RefDataBuffer_p[IndexRef4EC[m]].im(),
					-MicDataBuffer_p[IndexMic4EC[m]].re() * RefDataBuffer_p[IndexRef4EC[m]].im()
					+ MicDataBuffer_p[IndexMic4EC[m]].im() * RefDataBuffer_p[IndexRef4EC[m]].re());
				CorrTmpc[k] = VP_CPX(CorrTmpc[k].re() + tmp.re(), CorrTmpc[k].im() + tmp.im());
			}

			CorrTmp1[k] *= CorrTmp2[k];
			CorrTmp2[k] = CorrTmpc[k].re() * CorrTmpc[k].re() + CorrTmpc[k].im() * CorrTmpc[k].im();
		}

		for( int k = DetectionBins_e; k < FftBins_; ++k) {
			VP_CPX* MicDataBuffer_p = MicDataBuffer + k * MICNUM;
			VP_CPX* RefDataBuffer_p = RefDataBuffer + k * REFNUM;
			float* MicDataBuffer2_p = MicDataBufferAbs2 + k * MICNUM;
			float* RefDataBuffer2_p = RefDataBufferAbs2 + k * REFNUM;

			CorrTmp1[k] = MicDataBuffer2_p[IndexMic4EC[0]];
			CorrTmp2[k] = RefDataBuffer2_p[IndexRef4EC[0]];

			CorrTmpc[k] = VP_CPX(MicDataBuffer_p[IndexMic4EC[0]].re() * RefDataBuffer_p[IndexRef4EC[0]].re()
				+ MicDataBuffer_p[IndexMic4EC[0]].im() * RefDataBuffer_p[IndexRef4EC[0]].im(),
				-MicDataBuffer_p[IndexMic4EC[0]].re() * RefDataBuffer_p[IndexRef4EC[0]].im()
				+ MicDataBuffer_p[IndexMic4EC[0]].im() * RefDataBuffer_p[IndexRef4EC[0]].re());

			for (int m = 1; m < 3; ++m) {
				CorrTmp1[k] += MicDataBuffer2_p[IndexMic4EC[m]];
				CorrTmp2[k] += RefDataBuffer2_p[IndexRef4EC[m]];
				VP_CPX tmp = VP_CPX(MicDataBuffer_p[IndexMic4EC[m]].re() * RefDataBuffer_p[IndexRef4EC[m]].re()
					+ MicDataBuffer_p[IndexMic4EC[m]].im() * RefDataBuffer_p[IndexRef4EC[m]].im(),
					-MicDataBuffer_p[IndexMic4EC[m]].re() * RefDataBuffer_p[IndexRef4EC[m]].im()
					+ MicDataBuffer_p[IndexMic4EC[m]].im() * RefDataBuffer_p[IndexRef4EC[m]].re());
				CorrTmpc[k] = VP_CPX(CorrTmpc[k].re() + tmp.re(), CorrTmpc[k].im() + tmp.im());
			}

			CorrTmp1[k] *= CorrTmp2[k];
			CorrTmp2[k] = CorrTmpc[k].re() * CorrTmpc[k].re() + CorrTmpc[k].im() * CorrTmpc[k].im();
			SumUp += CorrTmp2[k];
			SumDn += CorrTmp1[k];
		}
		DetectionFlag[0] = CorrBand;
		for (int m = 1; m < FftBins_; ++m) {
			DetectionFlag[m] = (MicDataBufferAbs2[m * MICNUM + IndexMic4EC[0]] > RefDataBufferAbs2[m * REFNUM + IndexRef4EC[0]] * EC_BOUND_THRESH) ?
				0.0f : (CorrTmp2[m] / (CorrTmp1[m] + 1.0e-10f));
		}
	}
	else {
				for( int k = 1; k < DetectionBins_s; ++k) {
			VP_CPX* MicDataBuffer_p = MicDataBuffer + k * MICNUM;
			VP_CPX* RefDataBuffer_p = RefDataBuffer + k * REFNUM;
			float* MicDataBuffer2_p = MicDataBufferAbs2 + k * MICNUM;
			float* RefDataBuffer2_p = RefDataBufferAbs2 + k * REFNUM;

			CorrTmp1[k] = MicDataBuffer2_p[IndexMic4EC[0]];
			CorrTmp2[k] = RefDataBuffer2_p[IndexRef4EC[0]];

			CorrTmpc[k] = VP_CPX(MicDataBuffer_p[IndexMic4EC[0]].re() * RefDataBuffer_p[IndexRef4EC[0]].re()
				+ MicDataBuffer_p[IndexMic4EC[0]].im() * RefDataBuffer_p[IndexRef4EC[0]].im(),
				-MicDataBuffer_p[IndexMic4EC[0]].re() * RefDataBuffer_p[IndexRef4EC[0]].im()
				+ MicDataBuffer_p[IndexMic4EC[0]].im() * RefDataBuffer_p[IndexRef4EC[0]].re());

			for (int m = 1; m < 3; ++m) {
				CorrTmp1[k] += MicDataBuffer2_p[IndexMic4EC[m]];
				CorrTmp2[k] += RefDataBuffer2_p[IndexRef4EC[m]];
				VP_CPX tmp = VP_CPX(MicDataBuffer_p[IndexMic4EC[m]].re() * RefDataBuffer_p[IndexRef4EC[m]].re()
					+ MicDataBuffer_p[IndexMic4EC[m]].im() * RefDataBuffer_p[IndexRef4EC[m]].im(),
					-MicDataBuffer_p[IndexMic4EC[m]].re() * RefDataBuffer_p[IndexRef4EC[m]].im()
					+ MicDataBuffer_p[IndexMic4EC[m]].im() * RefDataBuffer_p[IndexRef4EC[m]].re());
				CorrTmpc[k] = VP_CPX(CorrTmpc[k].re() + tmp.re(), CorrTmpc[k].im() + tmp.im());
			}

			CorrTmp1[k] *= CorrTmp2[k];
			CorrTmp2[k] = CorrTmpc[k].re() * CorrTmpc[k].re() + CorrTmpc[k].im() * CorrTmpc[k].im();
		}

		for( int k = DetectionBins_e; k < FftBins_; ++k) {
			VP_CPX* MicDataBuffer_p = MicDataBuffer + k * MICNUM;
			VP_CPX* RefDataBuffer_p = RefDataBuffer + k * REFNUM;
			float* MicDataBuffer2_p = MicDataBufferAbs2 + k * MICNUM;
			float* RefDataBuffer2_p = RefDataBufferAbs2 + k * REFNUM;

			CorrTmp1[k] = MicDataBuffer2_p[IndexMic4EC[0]];
			CorrTmp2[k] = RefDataBuffer2_p[IndexRef4EC[0]];

			CorrTmpc[k] = VP_CPX(MicDataBuffer_p[IndexMic4EC[0]].re() * RefDataBuffer_p[IndexRef4EC[0]].re()
				+ MicDataBuffer_p[IndexMic4EC[0]].im() * RefDataBuffer_p[IndexRef4EC[0]].im(),
				-MicDataBuffer_p[IndexMic4EC[0]].re() * RefDataBuffer_p[IndexRef4EC[0]].im()
				+ MicDataBuffer_p[IndexMic4EC[0]].im() * RefDataBuffer_p[IndexRef4EC[0]].re());

			for (int m = 1; m < 3; ++m) {
				CorrTmp1[k] += MicDataBuffer2_p[IndexMic4EC[m]];
				CorrTmp2[k] += RefDataBuffer2_p[IndexRef4EC[m]];
				VP_CPX tmp = VP_CPX(MicDataBuffer_p[IndexMic4EC[m]].re() * RefDataBuffer_p[IndexRef4EC[m]].re()
					+ MicDataBuffer_p[IndexMic4EC[m]].im() * RefDataBuffer_p[IndexRef4EC[m]].im(),
					-MicDataBuffer_p[IndexMic4EC[m]].re() * RefDataBuffer_p[IndexRef4EC[m]].im()
					+ MicDataBuffer_p[IndexMic4EC[m]].im() * RefDataBuffer_p[IndexRef4EC[m]].re());
				CorrTmpc[k] = VP_CPX(CorrTmpc[k].re() + tmp.re(), CorrTmpc[k].im() + tmp.im());
			}

			CorrTmp1[k] *= CorrTmp2[k];
			CorrTmp2[k] = CorrTmpc[k].re() * CorrTmpc[k].re() + CorrTmpc[k].im() * CorrTmpc[k].im();
			SumUp += CorrTmp2[k];
			SumDn += CorrTmp1[k];
		}
		DetectionFlag[0] = CorrBand;
		for (int m = 1; m < FftBins_; ++m) {
			DetectionFlag[m] = (MicDataBufferAbs2[m * MICNUM + IndexMic4EC[0]] > RefDataBufferAbs2[m * REFNUM + IndexRef4EC[0]] * EC_BOUND_THRESH) ?
				0.0f : 0.83f * (CorrTmp2[m] / (CorrTmp1[m] + 1.0e-10f));
		}
	}
}

void VPAEC::EchoCancelingProcess(VP_CPX* OutData, VP_CPX* EchoEst, int IndexAlign)
{
	for (int m = 0; m < NumStage_a; ++m) {
		IndexRef4EC[m] = GetRefIndex(m + IndexAlign);
	}

	DetectModule();

	VP_CPX* Wa_ = Wa + IndexAlign * NumStage_a * FftBins_;
	VP_CPX* Pmr_ = Pmr + IndexAlign * NumStage_a * FftBins_;
	float* Prr_ = Prr + IndexAlign * NumStage_a * FftBins_;

	float beta = 0.07f;
	float beta_1 = 1.0f - beta;

	OutData[0] = 0.0f;
	for (int m = 1; m < IndexLowFreq; ++m) {
		VP_CPX* Wa_m = Wa_ + m * NumStage_a;
		VP_CPX* Pmr_m = Pmr_ + m * NumStage_a;
		float* Prr_m = Prr_ + m * NumStage_a;

		VP_CPX* MicDataBuffer_p = MicDataBuffer + m * MICNUM;
		VP_CPX* RefDataBuffer_p = RefDataBuffer + m * REFNUM;
		float* RefDataBuff2_p = RefDataBufferAbs2 + m * REFNUM;

		EchoEst[m] = 0.0f;
		if (DetectionFlag[m] > ThreshBins) {
			for (int k = 0; k < NumStage_a; ++k) {
				VP_CPX tmp1 = VP_CPX(MicDataBuffer_p[IndexMic4EC[0]].re() - EchoEst[m].re(),
					MicDataBuffer_p[IndexMic4EC[0]].im() - EchoEst[m].im());

				VP_CPX tpER = VP_CPX(tmp1.re() * RefDataBuffer_p[IndexRef4EC[k]].re()
					+ tmp1.im() * RefDataBuffer_p[IndexRef4EC[k]].im(),
					-tmp1.re() * RefDataBuffer_p[IndexRef4EC[k]].im()
					+ tmp1.im() * RefDataBuffer_p[IndexRef4EC[k]].re());

				float tpRR = RefDataBuff2_p[IndexRef4EC[k]];

				Pmr_m[k] = VP_CPX(Pmr_m[k].re() * beta_1, Pmr_m[k].im() * beta_1);
				Pmr_m[k] = VP_CPX(Pmr_m[k].re() + tpER.re() * beta, Pmr_m[k].im() + tpER.im() * beta);
				Prr_m[k] = Prr_m[k] * beta_1 + tpRR * beta;

				float tmp2 = Prr_m[k] + 1.0e-10f;
				Wa_m[k] = VP_CPX(Pmr_m[k].re() / tmp2, Pmr_m[k].im() / tmp2);

				EchoEst[m] = VP_CPX(EchoEst[m].re() + Wa_m[k].re() * RefDataBuffer_p[IndexRef4EC[k]].re()
					- Wa_m[k].im() * RefDataBuffer_p[IndexRef4EC[k]].im(),
					EchoEst[m].im() + Wa_m[k].re() * RefDataBuffer_p[IndexRef4EC[k]].im()
					+ Wa_m[k].im() * RefDataBuffer_p[IndexRef4EC[k]].re());
			}
		}
		else {
			for (int k = 0; k < NumStage_a; ++k) {
				EchoEst[m] = VP_CPX(EchoEst[m].re() + Wa_m[k].re() * RefDataBuffer_p[IndexRef4EC[k]].re()
					- Wa_m[k].im() * RefDataBuffer_p[IndexRef4EC[k]].im(),
					EchoEst[m].im() + Wa_m[k].re() * RefDataBuffer_p[IndexRef4EC[k]].im()
					+ Wa_m[k].im() * RefDataBuffer_p[IndexRef4EC[k]].re());	
			}
		}

		OutData[m] = VP_CPX(MicDataBuffer_p[IndexMic4EC[0]].re() - EchoEst[m].re(),
			MicDataBuffer_p[IndexMic4EC[0]].im() - EchoEst[m].im());

	}
 

	float beta = 0.15f;
	float beta_1 = 1.0f - beta;

	for (int m = IndexLowFreq; m < IndexMidFreq; ++m) {
		VP_CPX* Wa_m = Wa_ + m * NumStage_a;
		VP_CPX* Pmr_m = Pmr_ + m * NumStage_a;
		float* Prr_m = Prr_ + m * NumStage_a;

		VP_CPX* MicDataBuffer_p = MicDataBuffer + m * MICNUM;
		VP_CPX* RefDataBuffer_p = RefDataBuffer + m * REFNUM;
		float* RefDataBuff2_p = RefDataBufferAbs2 + m * REFNUM;

		EchoEst[m] = 0.0f;
		if (DetectionFlag[m] > ThreshBins) {
			for (int k = 0; k < NumStage_a; ++k) {
				VP_CPX tmp1 = VP_CPX(MicDataBuffer_p[IndexMic4EC[0]].re() - EchoEst[m].re(),
					MicDataBuffer_p[IndexMic4EC[0]].im() - EchoEst[m].im());

				VP_CPX tpER = VP_CPX(tmp1.re() * RefDataBuffer_p[IndexRef4EC[k]].re()
					+ tmp1.im() * RefDataBuffer_p[IndexRef4EC[k]].im(),
					-tmp1.re() * RefDataBuffer_p[IndexRef4EC[k]].im()
					+ tmp1.im() * RefDataBuffer_p[IndexRef4EC[k]].re());

				float tpRR = RefDataBuff2_p[IndexRef4EC[k]];

				Pmr_m[k] = VP_CPX(Pmr_m[k].re() * beta_1, Pmr_m[k].im() * beta_1);
				Pmr_m[k] = VP_CPX(Pmr_m[k].re() + tpER.re() * beta, Pmr_m[k].im() + tpER.im() * beta);
				Prr_m[k] = Prr_m[k] * beta_1 + tpRR * beta;

				float tmp2 = Prr_m[k] + 1.0e-10f;
				Wa_m[k] = VP_CPX(Pmr_m[k].re() / tmp2, Pmr_m[k].im() / tmp2);

				EchoEst[m] = VP_CPX(EchoEst[m].re() + Wa_m[k].re() * RefDataBuffer_p[IndexRef4EC[k]].re()
					- Wa_m[k].im() * RefDataBuffer_p[IndexRef4EC[k]].im(),
					EchoEst[m].im() + Wa_m[k].re() * RefDataBuffer_p[IndexRef4EC[k]].im()
					+ Wa_m[k].im() * RefDataBuffer_p[IndexRef4EC[k]].re());
			}
		}
		else {
			for (int k = 0; k < NumStage_a; ++k) {
				EchoEst[m] = VP_CPX(EchoEst[m].re() + Wa_m[k].re() * RefDataBuffer_p[IndexRef4EC[k]].re()
					- Wa_m[k].im() * RefDataBuffer_p[IndexRef4EC[k]].im(),
					EchoEst[m].im() + Wa_m[k].re() * RefDataBuffer_p[IndexRef4EC[k]].im()
					+ Wa_m[k].im() * RefDataBuffer_p[IndexRef4EC[k]].re());	
			}
		}

		OutData[m] = VP_CPX(MicDataBuffer_p[IndexMic4EC[0]].re() - EchoEst[m].re(),
			MicDataBuffer_p[IndexMic4EC[0]].im() - EchoEst[m].im());

	}

	float beta = 0.2f;
	float beta_1 = 1.0f - beta;

	for (int m = 1; m < IndexLowFreq; ++m) {
		VP_CPX* Wa_m = Wa_ + m * NumStage_a;
		VP_CPX* Pmr_m = Pmr_ + m * NumStage_a;
		float* Prr_m = Prr_ + m * NumStage_a;

		VP_CPX* MicDataBuffer_p = MicDataBuffer + m * MICNUM;
		VP_CPX* RefDataBuffer_p = RefDataBuffer + m * REFNUM;
		float* RefDataBuff2_p = RefDataBufferAbs2 + m * REFNUM;

		EchoEst[m] = 0.0f;
		if (DetectionFlag[m] > ThreshBins) {
			for (int k = 0; k < NumStage_a; ++k) {
				VP_CPX tmp1 = VP_CPX(MicDataBuffer_p[IndexMic4EC[0]].re() - EchoEst[m].re(),
					MicDataBuffer_p[IndexMic4EC[0]].im() - EchoEst[m].im());

				VP_CPX tpER = VP_CPX(tmp1.re() * RefDataBuffer_p[IndexRef4EC[k]].re()
					+ tmp1.im() * RefDataBuffer_p[IndexRef4EC[k]].im(),
					-tmp1.re() * RefDataBuffer_p[IndexRef4EC[k]].im()
					+ tmp1.im() * RefDataBuffer_p[IndexRef4EC[k]].re());

				float tpRR = RefDataBuff2_p[IndexRef4EC[k]];

				Pmr_m[k] = VP_CPX(Pmr_m[k].re() * beta_1, Pmr_m[k].im() * beta_1);
				Pmr_m[k] = VP_CPX(Pmr_m[k].re() + tpER.re() * beta, Pmr_m[k].im() + tpER.im() * beta);
				Prr_m[k] = Prr_m[k] * beta_1 + tpRR * beta;

				float tmp2 = Prr_m[k] + 1.0e-10f;
				Wa_m[k] = VP_CPX(Pmr_m[k].re() / tmp2, Pmr_m[k].im() / tmp2);

				EchoEst[m] = VP_CPX(EchoEst[m].re() + Wa_m[k].re() * RefDataBuffer_p[IndexRef4EC[k]].re()
					- Wa_m[k].im() * RefDataBuffer_p[IndexRef4EC[k]].im(),
					EchoEst[m].im() + Wa_m[k].re() * RefDataBuffer_p[IndexRef4EC[k]].im()
					+ Wa_m[k].im() * RefDataBuffer_p[IndexRef4EC[k]].re());
			}
		}
		else {
			for (int k = 0; k < NumStage_a; ++k) {
				EchoEst[m] = VP_CPX(EchoEst[m].re() + Wa_m[k].re() * RefDataBuffer_p[IndexRef4EC[k]].re()
					- Wa_m[k].im() * RefDataBuffer_p[IndexRef4EC[k]].im(),
					EchoEst[m].im() + Wa_m[k].re() * RefDataBuffer_p[IndexRef4EC[k]].im()
					+ Wa_m[k].im() * RefDataBuffer_p[IndexRef4EC[k]].re());	
			}
		}
        
		// get output of EC
		OutData[m] = VP_CPX(MicDataBuffer_p[IndexMic4EC[0]].re() - EchoEst[m].re(),
			MicDataBuffer_p[IndexMic4EC[0]].im() - EchoEst[m].im());
	}
}

void VPAEC::EchoCancelingProcess4Search()
{
	/* valid frame */
	{
		VectorClear(SearchFrameValid, SEARCHNUM);

		for (int m = 0; m < SEARCHNUM; ++m) {
			FeatureCopy[m] = FeatureCurrent[m];
		}

		int NumS = Min(EC_REF_ALIGN_MAX, 15);
		for (int m = 0; m < NumS; ++m) {
			int Idexmax = FindMax(FeatureCopy, SEARCHNUM);
			FeatureCopy[Idexmax] = -1.0f;
			SearchFrameValid[Idexmax] = 1;
		}

		for (int m = 0; m < SEARCHNUM; ++m) {
			FeatureCopy[m] = CancelOutSumSmth[m];
		}

		NumS = Min(EC_REF_ALIGN_MAX, 10);
		for (int m = 0; m < NumS; ++m) {
			int Indexmax = FindMax(FeatureCopy, SEARCHNUM);
			FeatureCopy[Indexmax] = -1.0f;
			SearchFrameValid[Indexmax] = 1;
		}

		for (int m = 0; m < SEARCHNUM; ++m) {
			FeatureCopy[m] = FeatureCopy[m];
		}

		for (int m = 0; m < 5; ++m) {
			int IdexMax = FindMax(FeatureCopy, SEARCHNUM);
			FeatureCopy[IdexMax] = -1.0f;
			SearchFrameValid[IdexMax] = 1;
		}
	}

	float beta = 0.05f;
	float beta_1 = 1.0f - beta;

	for (int frame = 0; frame < SEARCHNUM; ++frame) {
		CancelOutSumCurr[frame] = 0.0f;
		if (SearchFrameValid[frame] == 1) {
			for (int m = 0; m < 3; ++m) {
				IndexRef4EC[m] = GetRefIndex(m + frame);
			}

			DetectModule4Search();

			VP_CPX* Wa_ = Wa_s + frame * NumStage_s * SearchNum;
			VP_CPX* Pmr_ = Pmr_s + frame * NumStage_s * SearchNum;
			float*  Prr_ = Prr_s + frame * NumStage_s * SearchNum;

			for (int km = 0; km < SearchNum; ++km) {
				int k = km * 2 + SearchStart;
				VP_CPX* Wa_m = Wa_ + km * NumStage_s;

				VP_CPX* MicDataBuffer_p = MicDataBuffer + k * MICNUM;
				VP_CPX* RefDataBuffer_p = RefDataBuffer + k * REFNUM;
				float* MicDataBuffer2_p = MicDataBufferAbs2 + k * MICNUM;
				float* RefDataBuffer2_p = RefDataBufferAbs2 + k * REFNUM;

				EchoEst[k] = 0.0f;
				if (DetectionFlag[k] > ThreshBins) {
					VP_CPX* Pmr_m = Pmr_ + km * NumStage_s;
					float* Prr_m = Prr_ + km * NumStage_s;
					for (int j = 0; j < NumStage_s; ++j) {
						VP_CPX tmp1 = VP_CPX(MicDataBuffer_p[IndexMic4EC[0]].re() - EchoEst[k].re(),
							MicDataBuffer_p[IndexMic4EC[0]].im() - EchoEst[k].im());

						VP_CPX tpER = VP_CPX(tmp1.re() * RefDataBuffer_p[IndexRef4EC[j]].re()
							+ tmp1.im() * RefDataBuffer_p[IndexRef4EC[j]].im(),
							- tmp1.re() * RefDataBuffer_p[IndexRef4EC[j]].im()
							+ tmp1.im() * RefDataBuffer_p[IndexRef4EC[j]].re());

						float tpRR = RefDataBuffer2_p[IndexRef4EC[j]];

						Pmr_m[j] = VP_CPX(Pmr_m[j].re() * beta_1, Pmr_m[j].im() * beta_1);
						Pmr_m[j] = VP_CPX(Pmr_m[j].re() + tpER.re() * beta, Pmr_m[j].im() + tpER.im() * beta);
						Prr_m[j] = Prr_m[j] * beta_1 + tpRR * beta;

						float tmp2 = Prr_m[j] + 1.0e-10f;
						Wa_m[j] = VP_CPX(Pmr_m[j].re() / tmp2, Pmr_m[j].im() / tmp2);

						EchoEst[k] = VP_CPX(EchoEst[k].re() + Wa_m[j].re() * RefDataBuffer_p[IndexRef4EC[j]].re()
							- Wa_m[j].im() * RefDataBuffer_p[IndexRef4EC[j]].im(),
							EchoEst[k].im() + Wa_m[j].re() * RefDataBuffer_p[IndexRef4EC[j]].im()
							+ Wa_m[j].im() * RefDataBuffer_p[IndexRef4EC[j]].re());

					}
				}
				else {
					for (int j = 0; j < NumStage_s; ++j) {
						EchoEst[k] = VP_CPX(EchoEst[k].re() + Wa_m[j].re() * RefDataBuffer_p[IndexRef4EC[j]].re()
							- Wa_m[j].im() * RefDataBuffer_p[IndexRef4EC[j]].im(),
							EchoEst[k].im() + Wa_m[j].re() * RefDataBuffer_p[IndexRef4EC[j]].im()
							+ Wa_m[j].im() * RefDataBuffer_p[IndexRef4EC[j]].re());
					}
				}

				VP_CPX OutData = VP_CPX(MicDataBuffer_p[IndexMic4EC[0]].re() - EchoEst[k].re(),
					MicDataBuffer_p[IndexMic4EC[0]].im() - EchoEst[k].im());

				float cancel_out = MicDataBuffer2_p[IndexMic4EC[0]] - OutData.re() * OutData.re() - OutData.im() * OutData.im();
				cancel_out = (cancel_out > 0.0f) ? cancel_out : 0.0f;
				CancelOutSumCurr[frame] += cancel_out;

			}
		}
	}

	float OutSumSrc = 0.0f;
	for (int km = 0; km < SearchNum; ++km) {
		int k = km * 2 + SearchStart;
		OutSumSrc += MicDataBufferAbs2[k * MICNUM + IndexMic4EC[0]];
	}
	
	int IndexCancelMax = FindMax(CancelOutSumCurr, SEARCHNUM);

	if (CancelOutSumCurr[IndexCancelMax] > (0.2f * OutSumSrc + 1.0e-4f)) {
		for (int frame = 0; frame < SEARCHNUM; ++frame) {
			CancelOutSumCurr[frame] = CancelOutSumCurr[frame] * 0.8f + CancelOutSumCurr[frame] * 0.2f;
		}

		IndexAlignSearch = FindMax(CancelOutSumSmth, SEARCHNUM);
	}
}