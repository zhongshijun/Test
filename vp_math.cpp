#include "vp_math.h"

void VectorClear(int16_t *data, const int32_t num)
{
	for (int i = 0; i < num; ++i) {
		data[i] = 0;
	}
}

void VectorClear(int32_t* data, const int32_t num)
{
	for (int i = 0; i < num; ++i) {
		data[i] = 0;
	}
}

void VectorClear(float *data, const int32_t num)
{
	for (int i = 0; i < num; ++i) {
		data[i] = 0.f;
	}
}

void VectorClear(VP_CPX* data, const int32_t num)
{
	for (int i = 0; i < num; ++i) {
		data[i] = 0.f;
	}
}

void VectorClear(int32_t *data, const int32_t num)
{
	for (int i = 0; i < num; ++i) {
		data[i] = 0;
	}
}

int32_t FindMax(float *data, const int32_t num)
{
	int32_t Idx = 0;
	float dataMax = data[0];
	for (int i = 0; i < num; ++i) {
		if (dataMax < data[i]) {
			dataMax = data[i];
			Idx = i;
		}
	}
	return Idx;
}

void VectorAbs2(float *data2, const VP_CPX* data, const int32_t num)
{
	for (int i = 0; i < num; ++i) {
		data2[i] = data[i].re() * data[i].re() + data[i].im() * data[i].im();
	}
}






