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
	Float3& operator+=(const Float3& v) { x_ += v.x_; y_ += v.y_; z_ += v.z_; return *this; }
	Float3 operator*(const Float scalar) const { return Float3(x_ * scalar, y_ * scalar, z_ * scalar); }
	Float3& operator*=(const Float scalar) { x_ *= scalar; y_ *= scalar; z_ *= scalar; return *this; }
	Float Square() const { return x_ * x_ + y_ * y_ + z_ * z_; }
	Float Length() const { return std::sqrt(Square()); }
};
inline Float3 Normalize(const Float3& v)
{
	//Float length = v.Length();
	//Float inv_length = 1.0 / length;
	//return v * inv_length;
	Float d = std::sqrt(v.x_ * v.x_ + v.y_ * v.y_ + v.z_ * v.z_);
	return { v.x_ / d, v.y_ / d, v.z_ / d };
}
inline Float3 UniformSampleSphere(Float cosTheta, Float randForPhi)
{
	//return Normalize(Float3(cosTheta, randForPhi, z));
	auto sinTheta = std::sqrt(1 - cosTheta * cosTheta);
	auto phi = 2 * hx::PI * randForPhi;
	return Float3(sinTheta * std::cos(phi), sinTheta * std::sin(phi), cosTheta);
}
inline Float2 MapUpperSphereToPlane(const Float3& v)
{
	Float tmp = 1.0 / (std::abs(v.x_) + std::abs(v.y_) + std::abs(v.z_));
	return Float2(v.x_ * tmp, v.y_ * tmp);
}
inline Float3 MapPlaneToUpperSphere(const Float2& v)
{
	Float3 n;
	n.z_ = 1.0 - std::abs(v.x_) - std::abs(v.y_);
	n.x_ = v.x_;
	n.y_ = v.y_;
	n = Normalize(n);
	return n;
}
inline Float Ward(Float alpha, const Float3& wi, const Float3& wo)
{
	auto h = Normalize(wi + wo);
	auto cosThetaI = wi.z_;
	auto cosThetaO = wo.z_;
	auto tanThetaH = h.y_ / h.x_;
	auto tanThetaH2 = tanThetaH * tanThetaH;
	auto a2 = alpha * alpha;
	return std::exp(tanThetaH2 / a2) / (4.0 * PI * a2 * std::sqrt(cosThetaI * cosThetaO));
}
}