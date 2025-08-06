
#include "mathtypes.h"

namespace Ao {
	constexpr float DEFAULT_DEPTH = 64.0f;
	constexpr float DEFAULT_SCALE = 1.0f;
	constexpr float DEFAULT_GAMMA = 1.5f;

	void SetDepth(float depth = DEFAULT_DEPTH);
	void SetScale(float scale = DEFAULT_SCALE);
	void SetGamma(float gamma = DEFAULT_GAMMA);

	float Sample(float3_array pos, float3_array normal);
	float Blend(float src, float dest);
}; // namespace Ao
