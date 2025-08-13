
#include "mathtypes.h"

namespace Ao {
	constexpr float DEFAULT_DEPTH = 64.0f;
	constexpr float DEFAULT_INTENSITY = 1.0f;
	constexpr float DEFAULT_GAMMA = 1.5f;

	void SetDepth(float depth = DEFAULT_DEPTH);
	void SetIntensity(float scale = DEFAULT_INTENSITY);
	void SetGamma(float gamma = DEFAULT_GAMMA);

	float Sample(float3_array pos, float3_array normal);
	float Blend(float src, float dest);
}; // namespace Ao
