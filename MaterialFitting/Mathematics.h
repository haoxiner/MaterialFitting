#pragma once
#include <algorithm>
#include <numeric>
#include <cmath>

namespace hx
{
using Float = double;
constexpr Float PI = 3.14159265358979323846;
struct Float2
{
	Float x_, y_;
	Float2(Float scalar = 0.0) : Float2(scalar, scalar) {}
	Float2(Float x, Float y) : x_(x), y_(y) {}
};
struct Float3
{
	Float x_, y_, z_;
	Float3(Float scalar = 0.0) : Float3(scalar, scalar, scalar) {}
	Float3(Float x, Float y, Float z) : x_(x), y_(y), z_(z) {}
	Float3 operator+(const Float3& v) const { return Float3(x_ + v.x_, y_ + v.y_, z_ + v.z_); }
	Float3 operator*(const Float scalar) const { return Float3(x_ * scalar, y_ * scalar, z_ * scalar); }
	Float Square() const { return x_ * x_ + y_ * y_ + z_ * z_; }
	Float Length() const { return std::sqrt(Square()); }
};
inline Float3 Normalize(const Float3& v)
{
	Float length = v.Length();
	Float inv_length = 1.0 / length;
	return v * inv_length;
}
inline Float3 UniformSampleSphere(Float x, Float y)
{
	//return Normalize(Float3(x, y, z));
	auto sinTheta = std::sqrt(1 - x * x);
	auto phi = 2 * hx::PI * y;
	return Float3(sinTheta * std::cos(phi), sinTheta * std::sin(phi), x);
}
Float2 MapUpperSphereToPlane(const Float3& v)
{
	Float tmp = 1.0 / (std::abs(v.x_) + std::abs(v.y_) + std::abs(v.z_));
	return Float2(v.x_ * tmp, v.y_ * tmp);
}
}