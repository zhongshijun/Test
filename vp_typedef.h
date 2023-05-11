#include <stdint.h>

struct InterParamters
{
	float ProbFrame;

	bool ZeroFlag;
	bool EchoStateFlag;
	bool EchoStateFlagIm;

	float* PsigIn;
	float* Precho;
	float* PsigOut;
	float* Pnoise;
	float* Psig;

	float* GainEC;
	float* GainNS;
	float* GainAINR;
	float GainAGC;
};

struct Feature
{
	bool valid;
	uint32_t fe;
};