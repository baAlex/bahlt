/*

Copyright (c) 2025 Alexander Brandt

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the “Software”), to
deal in the Software without restriction, including without limitation the
rights to use, copy, modify, merge, publish, distribute, sublicense, and/or
sell copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED “AS IS”, WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS
IN THE SOFTWARE.

----

Or under public license. Code here is so simple that it will be a shame not
to share it properly. Put my name somewhere, if you want.

*/

#include "ao.hpp"

#include "hlrad.h"
#include "mathlib.h"

#include <cmath>

static float s_ao_depth = Ao::DEFAULT_DEPTH;
static float s_ao_scale = Ao::DEFAULT_SCALE;
static float s_ao_gamma = Ao::DEFAULT_GAMMA;

void Ao::SetDepth(float depth) {
	s_ao_depth = std::max(depth, 0.0f);
}

void Ao::SetScale(float scale) {
	s_ao_scale = std::max(scale, 0.0f);
}

void Ao::SetGamma(float gamma) {
	s_ao_gamma = std::max(gamma, 0.0f);
}

// Vertices of an icosphere hemisphere. Unlike other sphere algorithms its
// vertices are fairly equally distributed... at least that is what I known
// from Blender; regular sphere there have more vertices near the poles,
// which isn't desirable to use as a kernel
static constexpr int KERNEL_SAMPLES = 73;
static constexpr float3_array KERNEL[KERNEL_SAMPLES] = {
	{ 0, -0.525731086730957, 0.8506507873535156 },
	{ 0, 0.525731086730957, 0.8506507873535156 },
	{ 0.8506507873535156, 0, 0.525731086730957 },
	{ -0.8506507873535156, 0, 0.525731086730957 },
	{ -0.80901700258255, 0.5, 0.30901700258255005 },
	{ -0.5, 0.30901700258255005, 0.80901700258255 },
	{ -0.30901700258255005, 0.80901700258255, 0.5 },
	{ 0.30901700258255005, 0.80901700258255, 0.5 },
	{ -0.80901700258255, -0.5, 0.30901700258255005 },
	{ -0.5, -0.30901700258255005, 0.80901700258255 },
	{ 0, 0, 1 },
	{ 0.5, 0.30901700258255005, 0.80901700258255 },
	{ 0.80901700258255, 0.5, 0.30901700258255005 },
	{ 0.80901700258255, -0.5, 0.30901700258255005 },
	{ 0.5, -0.30901700258255005, 0.80901700258255 },
	{ 0.30901700258255005, -0.80901700258255, 0.5 },
	{ -0.30901700258255005, -0.80901700258255, 0.5 },
	{ -0.6817183494567871, 0.7165669202804565, 0.14762090146541595 },
	{ -0.5877852439880371, 0.6881909370422363, 0.4253253936767578 },
	{ -0.44286268949508667, 0.8641878366470337, 0.23885564506053925 },
	{ -0.7165669202804565, 0.14762090146541595, 0.6817183494567871 },
	{ -0.6881909370422363, 0.4253253936767578, 0.5877852439880371 },
	{ -0.8641878366470337, 0.23885564506053925, 0.44286268949508667 },
	{ -0.14762090146541595, 0.6817183494567871, 0.7165669202804565 },
	{ -0.4253253936767578, 0.5877852439880371, 0.6881909370422363 },
	{ -0.23885564506053925, 0.44286268949508667, 0.8641878366470337 },
	{ -0.1624598503112793, 0.9510565400123596, 0.2628655433654785 },
	{ 0.14762090146541595, 0.6817183494567871, 0.7165669202804565 },
	{ 0, 0.8506507873535156, 0.5257311463356018 },
	{ 0.1624598503112793, 0.9510565400123596, 0.2628655433654785 },
	{ 0.44286268949508667, 0.8641878366470337, 0.23885564506053925 },
	{ -0.9510565400123596, 0.2628655433654785, 0.1624598503112793 },
	{ -0.955422580242157, 0, 0.2952418029308319 },
	{ -0.9510565400123596, -0.2628655433654785, 0.1624598503112793 },
	{ -0.8641878366470337, -0.23885564506053925, 0.44286268949508667 },
	{ -0.6817183494567871, -0.7165669202804565, 0.14762090146541595 },
	{ -0.2628655433654785, 0.1624598503112793, 0.9510565400123596 },
	{ 0, 0.2952418029308319, 0.955422580242157 },
	{ -0.7165669202804565, -0.14762090146541595, 0.6817183494567871 },
	{ -0.5257311463356018, 0, 0.8506507873535156 },
	{ 0, -0.2952418029308319, 0.955422580242157 },
	{ -0.2628655433654785, -0.1624598503112793, 0.9510565400123596 },
	{ -0.23885564506053925, -0.44286268949508667, 0.8641878366470337 },
	{ 0.5877852439880371, 0.6881909370422363, 0.4253253936767578 },
	{ 0.6817183494567871, 0.7165669202804565, 0.14762090146541595 },
	{ 0.23885564506053925, 0.44286268949508667, 0.8641878366470337 },
	{ 0.4253253936767578, 0.5877852439880371, 0.6881909370422363 },
	{ 0.8641878366470337, 0.23885564506053925, 0.44286268949508667 },
	{ 0.6881909370422363, 0.4253253936767578, 0.5877852439880371 },
	{ 0.7165669202804565, 0.14762090146541595, 0.6817183494567871 },
	{ 0.6817183494567871, -0.7165669202804565, 0.14762090146541595 },
	{ 0.5877852439880371, -0.6881909370422363, 0.4253253936767578 },
	{ 0.44286268949508667, -0.8641878366470337, 0.23885564506053925 },
	{ 0.7165669202804565, -0.14762090146541595, 0.6817183494567871 },
	{ 0.6881909370422363, -0.4253253936767578, 0.5877852439880371 },
	{ 0.8641878366470337, -0.23885564506053925, 0.44286268949508667 },
	{ 0.14762090146541595, -0.6817183494567871, 0.7165669202804565 },
	{ 0.4253253936767578, -0.5877852439880371, 0.6881909370422363 },
	{ 0.23885564506053925, -0.44286268949508667, 0.8641878366470337 },
	{ 0.1624598503112793, -0.9510565400123596, 0.2628655433654785 },
	{ -0.14762090146541595, -0.6817183494567871, 0.7165669202804565 },
	{ 0, -0.8506507873535156, 0.5257311463356018 },
	{ -0.1624598503112793, -0.9510565400123596, 0.2628655433654785 },
	{ -0.44286268949508667, -0.8641878366470337, 0.23885564506053925 },
	{ 0.9510565400123596, -0.2628655433654785, 0.1624598503112793 },
	{ 0.955422580242157, 0, 0.2952418029308319 },
	{ 0.9510565400123596, 0.2628655433654785, 0.1624598503112793 },
	{ 0.2628655433654785, -0.1624598503112793, 0.9510565400123596 },
	{ 0.5257311463356018, 0, 0.8506507873535156 },
	{ 0.2628655433654785, 0.1624598503112793, 0.9510565400123596 },
	{ -0.5877852439880371, -0.6881909370422363, 0.4253253936767578 },
	{ -0.4253253936767578, -0.5877852439880371, 0.6881909370422363 },
	{ -0.6881909370422363, -0.4253253936767578, 0.5877852439880371 },
};

static float3_array sTBNRotation(float3_array vec, float3_array normal) {
	// Arbitrary vector that's not parallel to normal [a][d],
	// weird hack from a random fella on internet
	float3_array arbitrary = { 0.0f, 0.0f, 1.0f };

	if (std::abs(normal[1]) < 0.999f) {
		arbitrary[1] = 1.0f;
		arbitrary[2] = 0.0f;
	}

	// Calculate TBN [b][c], normal stuff
	float3_array tangent = cross_product(arbitrary, normal);
	normalize_vector(tangent);

	float3_array bi_tangent = cross_product(normal, tangent);

	// clang-format off
	float3_array ret;
	ret[0] = vec[0] * tangent[0] + vec[1] * bi_tangent[0] + vec[2] * normal[0];
	ret[1] = vec[0] * tangent[1] + vec[1] * bi_tangent[1] + vec[2] * normal[1];
	ret[2] = vec[0] * tangent[2] + vec[1] * bi_tangent[2] + vec[2] * normal[2];
	// clang-format on

	return ret;

	// Thanks random fella in a random forum:
	// [a] Erik Rufelt (2009): «There are infinitely many vectors that you
	// can use as tangents, as it can be any vector perpendicular to the
	// normal. You can find one by calculating the cross-product of the
	// normal and an arbitrary vector that isn't parallel to the normal.»
	// https://www.gamedev.net/forums/topic/552411-calculate-tangent-from-normal/

	// Here's something else with further Google-ing:
	// [d] Geeks3D (20013) Normal Mapping without Precomputed Tangent Space
	// Vectors
	// https://geeks3d.com/20130122/normal-mapping-without-precomputed-tangent-space-vectors

	// In any case is a weird hack because I'm not a graphical programmer

	// UPDATE, MORE SOURCES, ***IS INDEED A WEIRD HACK***:
	// Self Shadow, Perpendicular Possibilities (2011)
	// https://blog.selfshadow.com/2011/10/17/perp-vectors/
	// « A quick hack involves taking the cross product of the original unit
	// vector – let’s call it u(x,y,z) – with a fixed ‘up’ axis, e.g.
	// (0,1,0), and then normalising. A problem here is that if the two
	// vectors are very close – or equally, pointing directly away from each
	// other – then the result will be a degenerate vector. However, it’s
	// still a reasonable approach in the context of a camera, if the view
	// direction can be restricted to guard against this. A general solution
	// in this situation is to fall back to an alternative axis»

	// Is a rabbit hole going back to Hughes-Möller in the 90s, including
	// Pixar as well: Pixar, Building an Orthonormal Basis, Revisited
	// (2017). Journal of Computer Graphics Techniques Vol. 6, No. 1.

	// [b] https://shaderlabs.org/wiki/Shader_Tricks
	// [c] https://learnopengl.com/Advanced-Lighting/Normal-Mapping
}

float Ao::Sample(float3_array pos, float3_array normal) {
	float acc = 0.0f;

	for (int i = 0; i < KERNEL_SAMPLES; i += 1) {
		float3_array end_pos = sTBNRotation(KERNEL[i], normal);
		end_pos = vector_add(pos, vector_scale(end_pos, s_ao_depth));

		float frac;
		contents_t content = TestLineFrac(pos, end_pos, frac);

		if (content == contents_t::SOLID) {
			acc += frac;
		}
	}

	acc /= static_cast<float>(KERNEL_SAMPLES);
	acc = std::pow(acc, s_ao_gamma) * s_ao_scale;

	return 1.0f - std::clamp(acc, 0.0f, 1.0f);
}

float Ao::Blend(float src, float dest) {
	// Using a fancy blend function like 'Overlay' (Photoshop) looks like a
	// good idea, problem is that results are difficult to control and
	// inconsistent. That said, is a good idea to explore further.

	if (0) { // Developers, developers, developers
		return src * 128.0f;
	}

	return dest * src;
}
