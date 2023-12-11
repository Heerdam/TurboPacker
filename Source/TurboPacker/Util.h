#pragma once

#include "CoreMinimal.h"
#include "Kismet/GameplayStatics.h"
#include "DrawDebugHelpers.h"

#include <tuple>
#include <sstream>
#include <vector>
#include <array>
#include <array>
#include <queue>
#include <bitset>
#include <memory>
#include <cmath>
#include <optional>
#include <algorithm>
#include <type_traits>
#include <random>
#include <utility>
#include <type_traits>
#include <span>
#include <iostream>
#include <variant>
#include <complex>
#include <concepts>
#include <exception>
#include <chrono>
#include <ranges>
#include <filesystem>

#include <libmorton/morton.h>
#include <tsl/robin_map.h>
#include <tsl/robin_set.h>

#include <Eigen/Core>
#include <unsupported/Eigen/FFT>

#include "Logger.hpp"

using Dist = ::std::uniform_int_distribution<>;
using DistF = ::std::uniform_real_distribution<float>;
using DistD = ::std::uniform_real_distribution<double>;
using Rand = std::mt19937_64;

namespace Util {

	namespace Concepts {

		template <int32 DIM, typename T>
		concept IsEigenVector = requires(T v) {
			requires T::RowsAtCompileTime == DIM || T::ColsAtCompileTime == DIM;
		};

		template <typename T>
		concept IsUEVector2 = requires(T v) {
			requires std::is_same_v<T, UE::Math::TIntVector2<typename T::IntType>> ||
			std::is_same_v<T, UE::Math::TVector2<typename T::FReal>>;
		};

		template <typename T>
		concept IsUEVector3 = requires(T v) {
			requires std::is_same_v<T, UE::Math::TIntVector3<typename T::IntType>> ||
			std::is_same_v<T, UE::Math::TVector<typename T::FReal>>;
		};

		template<class TO, class FROM>
		[[nodiscard]] TO convert(const FROM& _vec) {
			static_assert((IsEigenVector<3, TO> || IsUEVector3<TO>) && (IsEigenVector<3, FROM> || IsUEVector3<FROM>));
			if constexpr (std::is_same_v<TO, FROM>) return _vec;
			else if constexpr (IsEigenVector<3, FROM>) return UE::Math::TVector<typename FROM::Scalar>(_vec(0), _vec(1), _vec(2));
			else return Eigen::Vector<typename FROM::FReal, 3>(_vec[0], _vec[1], _vec[2]);
		}//convert

		template<class V1, class V2>
		[[nodiscard]] float dist2(const V1& _v1, const V2& _v2) {
			static_assert((IsEigenVector<3, V1> || IsUEVector3<V1>) && (IsEigenVector<3, V2> || IsUEVector3<V2>));
			if constexpr (IsEigenVector<3, V1> && IsEigenVector<3, V2>)
				return (_v2 - _v1).squaredNorm();
			else if constexpr (IsUEVector3<V1> && IsUEVector3<V2>)
				return (_v2 - _v1).SquaredLength();
			else if constexpr (IsEigenVector<3, V1> && IsUEVector3<V2>)
				return FVector(_v2[0] - _v1(0), _v2[1] - _v1(1), _v2[2] - _v1(2)).SquaredLength();
			else return FVector(_v2(0) - _v1[0], _v2(1) - _v1[1], _v2(2) - _v1[2]).SquaredLength();
		}//dist2

		template <typename T>
		struct TypeIndex;

		template <>
		struct TypeIndex<float> {
			static const int value = 0;
		};

		template <>
		struct TypeIndex<double> {
			static const int value = 1;
		};

		template<class T>
		constexpr static int32 TypeIndex_v = TypeIndex<T>::value;

	}//Concepts

	//--------------------------------------

	template <class T>
	inline FVector c_to_HSL(const std::complex<T>& _c, T MAXVALUE = 10.) {
		const T H = std::clamp(std::abs(std::fmod(std::arg(_c), 2. * PI)), 0., 2. * PI);
		const T S = 1.;
		const T L = std::clamp(std::abs(MAXVALUE * std::atan(std::abs(_c)) / (0.5 * PI)), 0., 1.);
		return { H, S, L };
	}//c_to_HSL

	template <class T>
	inline FVector HSL_to_RGB_deg(const FVector& _hsl) {

		const T H = _hsl[0];
		const T S = _hsl[1];
		const T L = _hsl[2];

		assert(0. <= H && H <= 360.);
		assert(0. <= S && S <= 1.);
		assert(0. <= L && L <= 1.);

		const T C = (T(1.) - std::abs<T>(T(2.) * L - T(1.))) * S;
		const T X = C * (T(1.) - std::abs<T>(std::fmod<T>(H / T(60.), T(2.)) - T(1.)));
		const T m = L - C * T(0.5);

		//std::cout << H << ", " << S << ", " << L << std::endl;

		switch (size_t(H / 60.)) {
			case 0: return { C + m, X + m, m };
			case 1: return { X + m, C + m, m };
			case 2: return { m, C + m, X + m };
			case 3: return { m, X + m, C + m };
			case 4: return { X + m, m, C + m };
			case 5: return { C + m, m, X + m };
			default: return { 0., 0., 0. };
		}

	}//HSL_to_RGB_deg

	template <class T>
	inline FVector HSL_to_RGB_rad(const FVector& _hsl) {
		return HSL_to_RGB_deg<T>({ _hsl[0] * T(180. / PI), _hsl[1], _hsl[2]});
	};

	template<class T>
	[[nodiscard]] inline uint_fast64_t morton_index(T _x, T _y, T _z) {
		const uint_fast64_t out = libmorton::morton3D_64_encode(uint_fast64_t(std::nearbyint(std::abs(_x))), uint_fast64_t(std::nearbyint(std::abs(_y))), uint_fast64_t(std::nearbyint(std::abs(_z))));
		return (out & 0x1FFFFFFFFFFFFFFF) | uint_fast64_t(!std::signbit(double(_x))) << 63 | uint_fast64_t(!std::signbit(double(_y))) << 62 | uint_fast64_t(!std::signbit(double(_z))) << 61;
	}//morton_index_3d

	template<class T>
	[[nodiscard]] inline uint_fast64_t morton_index(T _x, T _y) {
		const uint_fast64_t out = libmorton::morton2D_64_encode(uint_fast64_t(std::nearbyint(std::abs(_x))), uint_fast64_t(std::nearbyint(std::abs(_y))));
		return (out & 0x1FFFFFFFFFFFFFFF) | uint_fast64_t(!std::signbit(double(_x))) << 63 | uint_fast64_t(!std::signbit(double(_y))) << 62;
	}//morton_index_2d

	template<class VECTOR>
		requires Concepts::IsEigenVector<3, VECTOR> || Concepts::IsUEVector3<VECTOR>
	[[nodiscard]] inline uint_fast64_t morton_index(const VECTOR& _x) {
		if constexpr (Concepts::IsEigenVector<3, VECTOR>)
			return morton_index(_x(0), _x(1), _x(2));
		else return morton_index(_x[0], _x[1], _x[2]);
	}//morton_index_3d

	template<class VECTOR>
		requires Concepts::IsEigenVector<2, VECTOR> || Concepts::IsUEVector2<VECTOR>
	[[nodiscard]] inline uint_fast64_t morton_index(const VECTOR& _x) {
		if constexpr (Concepts::IsEigenVector<2, VECTOR>)
			return morton_index(_x(0), _x(1));
		else return morton_index(_x[0], _x[1]);
	}//morton_index_3d

	////--------------------------------------

	template<class VECTOR>
	[[nodiscard]] inline VECTOR morton_index_2d(uint_fast64_t _idx) {
		const int32 sx = ((_idx >> 63) & 1) == 1 ? 1 : -1;
		const int32 sy = ((_idx >> 62) & 1) == 1 ? 1 : -1;
		uint_fast32_t x, y;
		libmorton::morton2D_64_decode(_idx, x, y);
		if constexpr (Concepts::IsEigenVector<2, VECTOR>) {
			using T = typename VECTOR::Scalar;
			return Eigen::Vector<typename VECTOR::Scalar, 2>(T(sx * int32(x)), T(sy * int32(y)));
		} else {
			using T = typename VECTOR::IntType;
			return UE::Math::TIntVector2<typename VECTOR::IntType>(T(sx * int32(x)), T(sy * int32(y)));
		}
	}//morton_index_3d

	template<class VECTOR>
	[[nodiscard]] inline VECTOR morton_index_3d(uint_fast64_t _idx) {
		const int32 sx = ((_idx >> 63) & 1) == 1 ? 1 : -1;
		const int32 sy = ((_idx >> 62) & 1) == 1 ? 1 : -1;
		const int32 sz = ((_idx >> 61) & 1) == 1 ? 1 : -1;
		uint_fast32_t x, y, z;
		libmorton::morton3D_64_decode(_idx, x, y, z);
		if constexpr (Concepts::IsEigenVector<3, VECTOR>) {
			using T = typename VECTOR::Scalar;
			return Eigen::Vector<typename VECTOR::Scalar, 3>(T(sx * int32(x)), T(sy * int32(y)), T(sz * int32(z)));
		} else {
			using T = typename VECTOR::IntType;
			return UE::Math::TIntVector3<typename VECTOR::IntType>(T(sx * int32(x)), T(sy * int32(y)), T(sz * int32(z)));
		}
	}//morton_index_3d

}//Util

namespace UE::Math {

	inline std::ostream& operator<< (std::ostream& s, const FVector& _r) {
		s << "[" << _r.X << ", " << _r.Y << ", " << _r.Z << "]";
		return s;
	}

	inline std::ostream& operator<< (std::ostream& s, const FIntVector& _r) {
		s << "[" << _r.X << ", " << _r.Y << ", " << _r.Z << "]";
		return s;
	}

	inline std::ostream& operator<< (std::ostream& s, const FInt32Point& _r) {
		s << "[" << _r.X << ", " << _r.Y << "]";
		return s;
	}

	inline std::ostream& operator<< (std::ostream& s, const FIntRect& _r) {
		s << _r.Min << _r.Max;
		return s;
	}

	inline std::ostream& operator<< (std::ostream& s, const FVector2D& _r) {
		s << "[" << _r.X << ", " << _r.Y << "]";
		return s;
	}

	inline std::ostream& operator<< (std::ostream& s, const FBox2D& _r) {
		s << _r.Min << _r.Max;
		return s;
	}

	inline std::ostream& operator<< (std::ostream& s, const FBox& _r) {
		s << _r.Min << _r.Max;
		return s;
	}

}
