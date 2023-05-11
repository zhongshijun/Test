#include <cstdint>
#include "vp_cpx.h"

#ifndef M_PI
#define M_PI 3.141592653589795f
#endif

#define B_01010101 (((uint32_t)(-1)) / 3)
#define B_00111001 (((uint32_t)(-1)) / 5)
#define B_00001111 (((uint32_t)(-1)) / 17)

void VectorClear(int16_t* data, const int32_t num);

void VectorClear(int32_t* data, const int32_t num);

void VectorClear(float* data, const int32_t num);

void VectorClear(VP_CPX* data, const int32_t num);

int32_t FindMax(float* data, const int32_t num);

void VectorAbs2(float *data2, const VP_CPX* data, const int32_t num);

inline float Min(float data1, float data2)
{
	return (data1 < data2) ? data1 : data2;
}

inline float Max(float data1, float data2)
{
	return (data1 > data2) ? data1 : data2;
}