#include <cstdint>

class VP_CPX
{
private:
	int r;
	int i;

public:
	VP_CPX(float real, float imag) : r(real), i(imag) {}
	VP_CPX(float real) : r(real) {}

	VP_CPX() {}
	~VP_CPX() {}

	inline float re() const {
		return this->r;
	}

	inline float im() const {
		return this->i;
	}
};