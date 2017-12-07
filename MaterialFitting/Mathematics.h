#pragma once
#include <algorithm>
#include <numeric>
#include <cmath>

namespace hx
{
using Float = double;
constexpr Float PI = 3.14159265358979323846;
Float DegreeToRadian(Float degree)
{
	return degree / 180.0 * PI;
}
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
	Float3& operator-(const Float3& v) const { return Float3(x_ - v.x_, y_ - v.y_, z_ - v.z_); }
	Float Square() const { return x_ * x_ + y_ * y_ + z_ * z_; }
	Float Length() const { return std::sqrt(Square()); }
	static Float Dot(const Float3& lhs, const Float3& rhs) { return lhs.x_ * rhs.x_ + lhs.y_ * rhs.y_ + lhs.z_ * rhs.z_; }
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
	//if (wi.z_ <= 0.0 || wo.z_ <= 0.0) {
	//	return 0.0;
	//}
	//auto h = Normalize(wi + wo);
	//auto cosThetaI = wi.z_;
	//auto cosThetaO = wo.z_;
	//
	//Float tanThetaH = 0.0;
	//
	//if (std::abs(h.x_) > 1e-10) {
	//	tanThetaH = h.y_ / h.x_;
	//}

	//auto tanThetaH2 = tanThetaH * tanThetaH;
	//auto a2 = alpha * alpha;
	//return std::exp(-tanThetaH2 / a2) / (4.0 * PI * a2 * std::sqrt(cosThetaI * cosThetaO));

	if (wi.z_ <= 0.0 || wo.z_ <= 0.0) {
		return 0.0;
	}
	auto H = Normalize(wi + wo);
	auto a2 = alpha * alpha;

	auto factor1 = 1.0 / (4.0 * PI * a2 * /*std::sqrt*/(wi.z_ * wo.z_));
	auto cosThetaH2 = H.z_ * H.z_;
	
	auto exponent = -(1.0 - cosThetaH2) / (cosThetaH2 * a2);

	Float specRef = factor1 * std::exp(exponent);

	//return specRef;
	// pdf * cosThetaI
	return specRef * wi.z_;
}
inline Float vMF(const Float3& mu, const Float k, const Float3& x)
{
	auto factor1 = k * Float3::Dot(mu, x);
	if (k < 100) {
		return k / 2.0 / PI / (std::exp(k - factor1) - std::exp(-k - factor1));
	} else {
		return k / 2.0 / PI * std::exp(k * (Float3::Dot(mu, x) - 1.0));
	}
}
inline Float3 SampleWard(Float alpha, Float u, Float v, const Float3& wo)
{
	Float phiH = std::atan(std::tan(2.0 * PI * v));
	//if (v > 0.5f)
	//	phiH += PI;
	Float cosPhiH = std::cos(phiH);
	Float sinPhiH = std::sqrt(1.0 - cosPhiH*cosPhiH);
	Float a2 = alpha * alpha;
	Float thetaH = std::atan(std::sqrt(-std::log(u) / (1.0 / (a2))));

	auto sinThetaH = std::sin(thetaH);
	auto cosThetaH = std::sqrt(1.0 - sinThetaH*sinThetaH);
	Float3 H = {sinThetaH * cosPhiH, sinThetaH * sinPhiH, cosThetaH};
	return H * (2.0 * Float3::Dot(wo, H)) - wo;
}
inline Float GGX(Float alpha, Float NdotH)
{
	auto a2 = alpha * alpha;
	auto NdotH2 = NdotH * NdotH;
	auto factor1 = (NdotH2 * (a2 - 1) + 1);
	if (factor1 < 1e-10) {
		return 1.0;
	}
	return a2 / (PI * factor1 * factor1);
}
}