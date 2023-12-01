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
#include <chrono>

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

	//--------------------------------------

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

	//--------------------------------------

	template<class T>
	class TURBOPACKER_API HeightMap {

		static_assert(sizeof(std::complex<float>) == sizeof(fftwf_complex));
		static_assert(sizeof(std::complex<double>) == sizeof(fftw_complex));
		static_assert(std::is_same_v<T, double> || std::is_same_v<T, float>);

	public:
		//real size
		const int32 x_;
		const int32 y_;
		const int32 z_;

		//padded size
		const int32 px_;
		const int32 py_;
		const int32 ls;

	private:
		//map
		std::vector<int32> map;

		//signal
		mutable std::vector<T> f_temp;
		mutable std::vector<std::complex<T>> c_temp;

		//kernel
		mutable std::vector<T> fi_temp;
		mutable std::vector<std::complex<T>> ci_temp;

		//temp
		mutable std::vector<std::complex<T>> si_temp;

		//result
		mutable std::array<std::vector<int32>, 6> result;

		//signal
		std::variant<fftwf_plan, fftw_plan> plan_r2c;
		std::variant<fftwf_plan, fftw_plan> plan_c2r;
		//kernel
		std::variant<fftwf_plan, fftw_plan> kplan_r2c;
		std::variant<fftwf_plan, fftw_plan> kplan_c2r;

	public:
		HeightMap(const int32 _x, const int32 _y, const int32 _z);
		~HeightMap();
		HeightMap(HeightMap&&) = default;
		HeightMap(const HeightMap&) = delete;
		//--------------------

		const std::array<std::vector<int32>, 6>&
			overlap(const FVector& _ext, const int32 _perm = 0x3F) const;

		void push(const FVector& _pos, const FVector& _ext);

		//--------------------
		// 
		//signal index
		int32 idx(const int32 _x, const int32 _y) const;

		//kernel index
		int32 kidx(const int32 _x, const int32 _y) const;

		//--------------------

		int32 operator[](const FIntVector2& _p) const { return map[idx(_p.X, _p.Y)]; }

		void print_size_in_bytes() const;

	};//HeightMap

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

//------------------------------------------------
//------------------ HeightMap -------------------
//------------------------------------------------

//pow(sqrt((w + 2 * (k/2 - 1)))

template<class T>
Util::HeightMap<T>::HeightMap(const int32 _x, const int32 _y, const int32 _z) : 
	x_(_x), y_(_y), z_(_z), 
	px_( std::pow(int32(std::sqrt(T(_x) + 2. * T(_x/2 - 1)) + 1), 2) ),
	py_( std::pow(int32(std::sqrt(T(_y) + 2. * T(_y/2 - 1)) + 1), 2) ),
	ls(px_ * py_)
{
	map.resize(_x * _y);

	//signal
	f_temp.resize(ls);
	c_temp.resize(ls);

	//kernel
	fi_temp.resize(_x * _y);
	ci_temp.resize(_x * _y);
	si_temp.resize(ls);

	//std::cout << px_ << "x" << py_ << std::endl;

	for(int32 i = 0; i < 6; ++i)
		result[i].resize(_x * _y);
	//-------------------------
	if constexpr (std::is_same_v<float, T>) {
		plan_r2c = rfftwf_plan_dft_r2c_2d(px_, py_, f_temp.data(), reinterpret_cast<fftwf_complex*>(c_temp.data()), FFTW_MEASURE);
		plan_c2r = fftwf_plan_dft_c2r_2d(px_, py_, reinterpret_cast<fftwf_complex*>(c_temp.data()), f_temp.data(), FFTW_MEASURE);
		kplan_r2c = fftwf_plan_dft_r2c_2d(x_, y_, fi_temp.data(), reinterpret_cast<fftwf_complex*>(ci_temp.data()), FFTW_MEASURE);
		kplan_c2r = fftwf_plan_dft_c2r_2d(x_, y_, reinterpret_cast<fftwf_complex*>(ci_temp.data()), fi_temp.data(), FFTW_MEASURE);
		ensure(std::get<0>(plan_r2c));
		ensure(std::get<0>(plan_c2r));
		ensure(std::get<0>(kplan_r2c));
		ensure(std::get<0>(kplan_c2r));
	} else {
		plan_r2c = fftw_plan_dft_r2c_2d(px_, py_, f_temp.data(), reinterpret_cast<fftw_complex*>(c_temp.data()), FFTW_MEASURE);
		plan_c2r = fftw_plan_dft_c2r_2d(px_, py_, reinterpret_cast<fftw_complex*>(c_temp.data()), f_temp.data(), FFTW_MEASURE);
		kplan_r2c = fftw_plan_dft_r2c_2d(x_, y_, fi_temp.data(), reinterpret_cast<fftw_complex*>(ci_temp.data()), FFTW_MEASURE);
		kplan_c2r = fftw_plan_dft_c2r_2d(x_, y_, reinterpret_cast<fftw_complex*>(ci_temp.data()), fi_temp.data(), FFTW_MEASURE);
		ensure(std::get<1>(plan_r2c));
		ensure(std::get<1>(plan_c2r));
		ensure(std::get<1>(kplan_r2c));
		ensure(std::get<1>(kplan_c2r));
	}
}//Util::HeightMap::HeightMap

template<class T>
Util::HeightMap<T>::~HeightMap() {
	if constexpr (std::is_same_v<float, T>) {
		fftwf_destroy_plan(std::get<0>(plan_r2c));
		fftwf_destroy_plan(std::get<0>(plan_c2r));
		fftwf_destroy_plan(std::get<0>(kplan_r2c));
		fftwf_destroy_plan(std::get<0>(kplan_c2r));
	} else {
		fftw_destroy_plan(std::get<1>(plan_r2c));
		fftw_destroy_plan(std::get<1>(plan_c2r));
		fftw_destroy_plan(std::get<1>(kplan_r2c));
		fftw_destroy_plan(std::get<1>(kplan_c2r));
	}
}//Util::HeightMap::~HeightMap

template<class T>
const std::array<std::vector<int32>, 6>&
Util::HeightMap<T>::overlap(const FVector& _ext, const int32 _perm) const {

	const T SCALE = T(1.) / T(ls);

	//fft map
	std::fill(f_temp.begin(), f_temp.end(), 0);
	for (int32 x = 0; x < x_; ++x) {
		for (int32 y = 0; y < y_; ++y) {
			const int32 i1 = kidx(x, y);
			const int32 i2 = idx(x, y);
			f_temp[i2] = T(map[i1]);
		}
	}

	//signal
	if constexpr (std::is_same_v<float, T>)
		fftwf_execute_dft_r2c(std::get<0>(plan_r2c), f_temp.data(), reinterpret_cast<fftwf_complex*>(c_temp.data()));
	else fftw_execute_dft_r2c(std::get<1>(plan_r2c), f_temp.data(), reinterpret_cast<fftw_complex*>(c_temp.data()));

	//permutations
	if (_perm | 0x1) {
		std::fill(fi_temp.begin(), fi_temp.end(), 0.f);

		const int32 w = int32(std::round(_ext[0])) / 2;
		const int32 h = int32(std::round(_ext[1])) / 2;
		const int32 cw = int32(x_ / 2);
		const int32 ch = int32(y_ / 2);

		//center kernel
		for (int32 x = cw - w; x <= cw + w; ++x) {
			for (int32 y = ch - h; y <= ch + h; ++y) {	
				const int32 i = kidx(x, y);
				fi_temp[i] = 1.f;
			}
		}
		/*
		for (int32 y = 0; y <= std::round(_ext[1]); ++y) {
			for (int32 x = 0; x <= std::round(_ext[0]); ++x) {
				const int32 i = idx(x, y);
				if (i < 0 || i >= map.size()) continue;
				fi_temp[i] = 1.f;
			}
		}
		*/

		if constexpr (std::is_same_v<float, T>)
			fftwf_execute_dft_r2c(std::get<0>(kplan_r2c), fi_temp.data(), reinterpret_cast<fftwf_complex*>(ci_temp.data()));
		else fftw_execute_dft_r2c(std::get<1>(kplan_r2c), fi_temp.data(), reinterpret_cast<fftw_complex*>(ci_temp.data()));
		
		for (int32 x = 0; x < x_; ++x) {
			for (int32 y = 0; y < y_; ++y) {
				const int32 i1 = kidx(x, y);
				const int32 i2 = idx(x, y);
				si_temp[i2] = c_temp[i2] * ci_temp[i1];
			}
		}

		if constexpr (std::is_same_v<float, T>)
			fftwf_execute_dft_c2r(std::get<0>(plan_c2r), reinterpret_cast<fftwf_complex*>(si_temp.data()), fi_temp.data());
		else fftw_execute_dft_c2r(std::get<1>(plan_c2r), reinterpret_cast<fftw_complex*>(si_temp.data()), fi_temp.data());

		std::fill(result[0].begin(), result[0].end(), 0);

		for (int32 x = 0; x < x_; ++x) {
			for (int32 y = 0; y < y_; ++y) {
				const int32 i1 = kidx(x, y);
				const int32 i2 = idx(x, y);

				const int32 expec = int32(std::round(_ext[0]) + 1) * int32(std::round(_ext[0]) + 1) * map[i1];
				const int32 r = int32(std::abs(std::round(fi_temp[i2] * SCALE)));
				result[0][i1] = r > expec ? -1 : r < map[i1] ? 1 : 0;
			}
		}

	}

	/*if (_perm | 0x2) {
	}
	
	if (_perm | 0x8) {
	}

	if (_perm | 0x10) {
	}

	if (_perm | 0x20) {
	}

	if (_perm | 0x40) {
	}*/

	return result;
	
}//Util::HeightMap::overlap

template<class T>
void Util::HeightMap<T>::push(const FVector& _pos, const FVector& _ext) {

	const int32 minX = int32(_pos[0] - _ext[0]);
	const int32 maxX = int32(_pos[0] + _ext[0]);

	const int32 minY = int32(_pos[1] - _ext[1]);
	const int32 maxY = int32(_pos[1] + _ext[1]);

	const int32 maxZ = int32(_pos[2] + _ext[2]);

	for (int32 x = minX; x <= maxX; ++x) {
		for (int32 y = minY; y <= maxY; ++y) {	
			const int32 i = idx(x, y);
			if (i < 0 || i >= map.size()) continue;//TODO height
			map[i] = maxZ;
		}
	}

}//Util::HeightMap::push

template<class T>
void Util::HeightMap<T>::print_size_in_bytes() const {
	std::cout << "Total size: " << double(map.size() * sizeof(int32) + f_temp.size() * sizeof(T) + c_temp.size() * sizeof(std::complex<T>) +
		fi_temp.size() * sizeof(T) + ci_temp.size() * sizeof(T) + 6 * result[0].size() * sizeof(int32)) / 1000000. << "Mb" << std::endl;
}//Util::HeightMap::print_size_in_bytes

template<class T>
int32 Util::HeightMap<T>::idx(const int32 _x, const int32 _y) const {
	const int32 ox = (px_ - x_) / 2;
	const int32 oy = (py_ - y_) / 2;
	//const int32 out = (_x + ox) + (_y + oy) * px_;
	const int32 out = (_y + oy) + (_x + ox) * py_;
	if (out < 0 || out >= f_temp.size()) {
		std::cout << "idx: " << out << std::endl;
		return 0;
	}
	return out;
}//Util::HeightMap::idx

template<class T>
int32 Util::HeightMap<T>::kidx(const int32 _x, const int32 _y) const {
	const int32 out = (_y) + (_x) * y_;
	if (out < 0 || out >= map.size()) {
		std::cout << "kidx: " << out << std::endl;
		return 0;
	}
	return out;
}//Util::HeightMap::idx
