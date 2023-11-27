#pragma once

#include "CoreMinimal.h"

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

namespace Util{

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

	}//Concepts

	//--------------------------------------

	template<class T>
	[[nodiscard]] uint_fast64_t morton_index(T _x, T _y, T _z) {
		const uint_fast64_t out = libmorton::morton3D_64_encode(uint_fast64_t(std::nearbyint(std::abs(_x))), uint_fast64_t(std::nearbyint(std::abs(_y))), uint_fast64_t(std::nearbyint(std::abs(_z))));
		return (out & 0x1FFFFFFFFFFFFFFF) | uint_fast64_t(!std::signbit(double(_x))) << 63 | uint_fast64_t(!std::signbit(double(_y))) << 62 | uint_fast64_t(!std::signbit(double(_z))) << 61;
	}//morton_index_3d

	template<class T>
	[[nodiscard]] uint_fast64_t morton_index(T _x, T _y) {
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

	//--------------------------------------

	template<class T>
	class TURBOPACKER_API InfGrid {

		static_assert(std::is_default_constructible_v<T>, "T is not default constructible");

		tsl::robin_map<uint_fast64_t, T> map;

		float dx, dy, dz;

	public:

		InfGrid() = default;

		template<class VECTOR3>
		InfGrid(const VECTOR3& _cell_size);

		template<class VECTOR>
		[[nodiscard]] T* init_cell(const VECTOR& _pos);

		template<class VECTOR>
		[[nodiscard]] T* operator[] (const VECTOR& _pos) const;

		[[nodiscard]] T* operator[] (uint_fast64_t _idx) const;

		[[nodiscard]] float get_dx() const { return dx; }
		[[nodiscard]] float get_dy() const { return dy; }
		[[nodiscard]] float get_dz() const { return dz; }

	};//InfGrid

	template<class T, class VECTOR>
	class TURBOPACKER_API InfBucketGrid {

	public:
		using VectorType = VECTOR;
		using PayloadType = T;
		using ReturnType = std::pair<T, VECTOR>;

	private:

		InfGrid<std::vector<ReturnType>> grid;

	public:

		InfBucketGrid() = default;

		template<class VECTOR3>
		InfBucketGrid(const VECTOR3& _cell_size);

		//-----------------
		template<class VEC>
		void push(const T& _val, const VEC& _point);
		//-----------------
		template<class VECTOR3>
		[[nodiscard]] std::vector<ReturnType*> radial_search(const VECTOR3& _pos, float _radius) const;
		//-----------------
		template<class VECTOR3>
		[[nodiscard]] ReturnType* nn(const VECTOR3& _point, float _max_search_distance = -1.f) const;
	};

}//Util

//------------------------------------------------
//-------------------InfGrid----------------------
//------------------------------------------------

template<class T>
template<class VECTOR3>
inline Util::InfGrid<T>::InfGrid(const VECTOR3& _cell_size) {
	if constexpr (Concepts::IsEigenVector<3, VECTOR3>) {
		dx = 1.f / float(_cell_size(0));
		dy = 1.f / float(_cell_size(1));
		dz = 1.f / float(_cell_size(2));
	} else {
		dx = 1.f / float(_cell_size[0]);
		dy = 1.f / float(_cell_size[1]);
		dz = 1.f / float(_cell_size[2]);
	}
}//Util::InfGrid::InfGrid

template<class T>
template<class VECTOR>
T* Util::InfGrid<T>::init_cell(const VECTOR& _pos) {
	if constexpr (Concepts::IsEigenVector<3, VECTOR>) {
		const uint_fast64_t idx = morton_index(VECTOR(_pos(0) * dx, _pos(1) * dy, _pos(2) * dz));
		map.emplace(idx, T());
		return &map[idx];
	} else {
		const uint_fast64_t idx = morton_index(VECTOR(_pos[0] * dx, _pos[1] * dy, _pos[2] * dz));
		map.emplace(idx, T());
		return &map[idx];
	}
}//Util::InfGrid::init_cell

template<class T>
template<class VECTOR>
inline T* Util::InfGrid<T>::operator[](const VECTOR& _pos) const {
	if constexpr (Concepts::IsEigenVector<3, VECTOR>) {
		const uint_fast64_t idx = morton_index(VECTOR(_pos(0) * dx, _pos(1) * dy, _pos(2) * dz));
		if (map.contains(idx)) return const_cast<T*>(&map.find(idx).value());
		return nullptr;
	} else {
		const uint_fast64_t idx = morton_index(VECTOR(_pos[0] * dx, _pos[1] * dy, _pos[2] * dz));
		if (map.contains(idx)) return const_cast<T*>(&map.find(idx).value());
		return nullptr;
	}
}//Util::InfGrid::operator[]

template<class T>
T* Util::InfGrid<T>::operator[] (uint_fast64_t _idx) const {
	if(map.contains(_idx)) return const_cast<T*>(&map.find(_idx).value());
	return nullptr;	
}//Util::InfGrid::operator[]

//------------------------------------------------
//-------------------InfBucketGrid----------------
//------------------------------------------------

template<class T, class VECTOR>
template<class VECTOR3>
inline Util::InfBucketGrid<T, VECTOR>::InfBucketGrid(const VECTOR3& _cell_size) {
	grid = InfGrid<std::vector<ReturnType>>(_cell_size);
}//Util::InfBucketGrid::InfBucketGrid

template<class T, class VECTOR>
template<class VEC>
inline void Util::InfBucketGrid<T, VECTOR>::push(const T& _val, const VEC& _point) {
	auto e = grid[_point];
	if (!e) e = grid.init_cell(_point);
	e->push_back({ _val, Concepts::convert<VECTOR, VEC>(_point) });
}//Util::InfBucketGrid::push

template<class T, class VECTOR>
template<class VECTOR3>
inline std::vector<typename Util::InfBucketGrid<T, VECTOR>::ReturnType*> 
Util::InfBucketGrid<T, VECTOR>::radial_search(const VECTOR3& _pos, float _radius) const {

	const int32 minX = int32((_pos[0] - _radius) * grid.get_dx());
	const int32 maxX = int32((_pos[0] + _radius) * grid.get_dx());

	const int32 minY = int32((_pos[1] - _radius) * grid.get_dy());
	const int32 maxY = int32((_pos[1] + _radius) * grid.get_dy());

	const int32 minZ = int32((_pos[2] - _radius) * grid.get_dz());
	const int32 maxZ = int32((_pos[2] + _radius) * grid.get_dz());

	const float r2 = _radius * _radius;

	std::vector<ReturnType*> out;
	for (int32 z = minZ; z <= maxZ; ++z) {
		for (int32 y = minY; y <= maxY; ++y) {
			for (int32 x = minX; x <= maxX; ++x) {

				const uint_fast64_t idx = morton_index(FVector(x, y, z));
				std::vector<ReturnType>* t = grid[idx];
				if (t) {
					for (auto& tt : (*t)) {
						const auto& [ty, p] = tt;
						const float d = Concepts::dist2(p, _pos);
						if (d <= _radius)
							out.push_back(&tt);
					}
				}

			}
		}
	}
	return out;
}//Util::InfBucketGrid::radial_search

template<class T, class VECTOR>
template<class VECTOR3>
inline typename Util::InfBucketGrid<T, VECTOR>::ReturnType* 
Util::InfBucketGrid<T, VECTOR>::nn(const VECTOR3& _point, float _max_search_distance) const {

	const float sd = _max_search_distance < 0.f ?
		2*std::max(std::max(grid.get_dx(), grid.get_dy()), grid.get_dz()) :
		_max_search_distance;

	const auto rs = radial_search(_point, sd);
	if (rs.empty()) return nullptr;
	ReturnType* out = nullptr;
	float dist = std::numeric_limits<float>::infinity();
	for (const auto& tt : rs) {
		const auto& [ty, p] = *tt;
		const float d = Concepts::dist2(p, _point);
		if (d < dist) {
			out = &(*tt);
			dist = d;
		}
	}

	return out;
}//Util::InfBucketGrid::nn
