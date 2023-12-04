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
#include <ranges>

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
		const int32 n0_;
		const int32 n1_;
		const int32 n1_c_;
		const int32 h_;

		//padded size
		const int32 p_n0_;
		const int32 p_n1_;
		const int32 p_n1_c_;

	private:
		//map
		std::vector<int32> map;

		//signal
		double* signal_r = nullptr;
		fftw_complex* signal_c = nullptr;

		//kernel
		double* kernel_r = nullptr;
		fftw_complex* kernel_c = nullptr;

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
		HeightMap(const int32 _n0, const int32 _n1, const int32 _y);
		~HeightMap();
		HeightMap(HeightMap&&) = default;
		HeightMap(const HeightMap&) = delete;
		//--------------------

		auto overlap(const FVector& _ext, const int32 _perm = 0x3F) const;

		void push(const FVector& _pos, const FVector& _ext);

		//--------------------
		// 
		//signal index
		int32 idx(const int32 _n0, const int32 _n1) const;

		//kernel index
		int32 kidx(const int32 _n0, const int32 _n1) const;

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
Util::HeightMap<T>::HeightMap(const int32 _n0, const int32 _n1, const int32 _h) : 
	n0_(_n0), n1_(_n1), n1_c_(_n1 / 2 + 1), h_(_h),
	p_n0_(2 * _n0 - 1), p_n1_(2 * _n1 - 1), p_n1_c_(p_n1_ / 2 + 1)
{
	map.resize(_n0 * _n1);
	std::fill(map.begin(), map.end(), 0);

	//signal
	signal_r = fftw_alloc_real(p_n0_ * p_n1_);
	signal_c = fftw_alloc_complex(p_n0_ * p_n1_c_);

	//kernel
	kernel_r = fftw_alloc_real(p_n0_ * p_n1_);
	kernel_c = fftw_alloc_complex(p_n0_ * p_n1_c_);

	//std::cout << px_ << "x" << py_ << std::endl;

	//for(int32 i = 0; i < 6; ++i)
		//result[i].resize(_x * _y);
	//-------------------------
	/*if constexpr (std::is_same_v<float, T>) {
		plan_r2c = fftwf_plan_dft_r2c_2d(py_, px_, f_temp.data(), reinterpret_cast<fftwf_complex*>(c_temp.data()), FFTW_MEASURE);
		plan_c2r = fftwf_plan_dft_c2r_2d(py_, px_, reinterpret_cast<fftwf_complex*>(c_temp.data()), f_temp.data(), FFTW_MEASURE);
		kplan_r2c = fftwf_plan_dft_r2c_2d(y_, x_, fi_temp.data(), reinterpret_cast<fftwf_complex*>(ci_temp.data()), FFTW_MEASURE);
		kplan_c2r = fftwf_plan_dft_c2r_2d(y_, x_, reinterpret_cast<fftwf_complex*>(ci_temp.data()), fi_temp.data(), FFTW_MEASURE);
		ensure(std::get<0>(plan_r2c));
		ensure(std::get<0>(plan_c2r));
		ensure(std::get<0>(kplan_r2c));
		ensure(std::get<0>(kplan_c2r));
	} else */
	{
		plan_r2c = fftw_plan_dft_r2c_2d(p_n0_, p_n1_, signal_r, signal_c, FFTW_MEASURE);
		plan_c2r = fftw_plan_dft_c2r_2d(p_n0_, p_n1_, signal_c, signal_r, FFTW_MEASURE);
		kplan_r2c = fftw_plan_dft_r2c_2d(p_n0_, p_n1_, kernel_r, kernel_c, FFTW_MEASURE);
		kplan_c2r = fftw_plan_dft_c2r_2d(p_n0_, p_n1_, kernel_c, kernel_r, FFTW_MEASURE);
		ensure(std::get<1>(plan_r2c));
		ensure(std::get<1>(plan_c2r));
		ensure(std::get<1>(kplan_r2c));
		ensure(std::get<1>(kplan_c2r));
	}
}//Util::HeightMap::HeightMap

template<class T>
Util::HeightMap<T>::~HeightMap() {
	/*if constexpr (std::is_same_v<float, T>) {
		fftwf_destroy_plan(std::get<0>(plan_r2c));
		fftwf_destroy_plan(std::get<0>(plan_c2r));
		fftwf_destroy_plan(std::get<0>(kplan_r2c));
		fftwf_destroy_plan(std::get<0>(kplan_c2r));
	} else*/ 
	{
		fftw_destroy_plan(std::get<1>(plan_r2c));
		fftw_destroy_plan(std::get<1>(plan_c2r));
		fftw_destroy_plan(std::get<1>(kplan_r2c));
		fftw_destroy_plan(std::get<1>(kplan_c2r));
		fftw_free(signal_r);
		fftw_free(signal_c);
		fftw_free(kernel_r);
		fftw_free(kernel_c);
	}
}//Util::HeightMap::~HeightMap

template<class T>
auto Util::HeightMap<T>::overlap(const FVector& _ext, const int32 _perm) const {

	/*
	 
		wp = (sqrt(w + (w + 1)) + 1)^2
		hp = (sqrt(h + (h + 1)) + 1)^2

		signal:
			in: wp x hp [n_0 x n_1]
			out: wp x (hp/2 + 1) [n_0 x (n_1/2 + 1)]

		kernel:
			in: w x h
			out: w x (h/2 + 1)

		convolution:
			in: wp x hp [n_0 x n_1]


	*/

	const int32 o1 = n1_ / 2;
	const int32 o1_c = n1_c_ / 2;
	const int32 o0 = n0_ / 2;

	const T SCALE = T(1.) / T(n0_ * n1_);

	for (int32 i = 0; i < p_n0_ * p_n1_; ++i)
		signal_r[i] = 0.;

	//fft map
	for (int32 n1 = 0; n1 < n1_; ++n1) {
		for (int32 n0 = 0; n0 < n0_; ++n0) {
			const int32 i1 = n1 + n0 * n1_;
			const int32 i2 = (n1 + o1) + (n0 + o0) * p_n1_;
			signal_r[i2] = T(map[i1]);
			//result[1][i1] = map[i1];
		}
	}

	fftw_execute(std::get<1>(plan_r2c));

	for (int32 i = 0; i < p_n0_ * p_n1_; ++i)
		kernel_r[i] = 0.;

	const int32 w = int32(std::round(_ext[0]));
	const int32 h = int32(std::round(_ext[1]));
	const int32 cw = int32(p_n0_ / 2);
	const int32 ch = int32(p_n1_ / 2);

	//center kernel
	/*for (int32 n1 = ch - h; n1 <= ch + h; ++n1) {	
		for (int32 n0 = cw - w; n0 <= cw + w; ++n0) {
			const int32 i = n1 + n0 * p_n1_;
			kernel_r[i] = 1.f;
		}
	}*/

	for (int32 n0 = 0; n0 <= h; ++n0) {
		for (int32 n1 = 0; n1 <= w; ++n1) {
			const int32 i = n1 + n0 * p_n1_;
			kernel_r[i] = 1.f;
		}
	}

	fftw_execute(std::get<1>(kplan_r2c));

	if constexpr (true) {
		for (int32 n0 = 0; n0 < p_n0_; ++n0) {
			for (int32 n1 = 0; n1 < p_n1_c_; ++n1) {
				const int32 i = n1 + n0 * p_n1_c_;
				const std::complex<T> c1(signal_c[i][0], signal_c[i][1]);
				const std::complex<T> c2(kernel_c[i][0], kernel_c[i][1]);
				const auto r = c1 * c2;
				signal_c[i][0] = r.real();
				signal_c[i][1] = r.imag();
			}
		}
	}

	fftw_execute(std::get<1>(plan_c2r));

	return signal_r;















	//signal

//if constexpr (std::is_same_v<float, T>)
//	fftwf_execute_dft_r2c(std::get<0>(plan_r2c), f_temp.data(), reinterpret_cast<fftwf_complex*>(c_temp.data()));
//else fftw_execute_dft_r2c(std::get<1>(plan_r2c), f_temp.data(), reinterpret_cast<fftw_complex*>(c_temp.data()));



	/*for (int32 x = 0; x < x_; ++x) {
		for (int32 y = 0; y < y_; ++y) {
			const int32 i1 = idx(x, y);
			const int32 i2 = kidx(x, y);
			result[1][i2] = int32(f_temp[i1]);
			result[1][i2] = int32(c_temp[i1].real());
		}
	}*/

	//fftw_execute(std::get<1>(plan_c2r));
	//for (int32 i = 0; i < f_temp.size(); ++i) {
	//	c_temp[i] = { f_temp[i], 0. };
		//std::cout << f_temp[i] << std::endl;
	//}

	

	////permutations
	//if (_perm | 0x1) {
	//	std::fill(fi_temp.begin(), fi_temp.end(), 0.f);

	//	const int32 w = int32(std::round(_ext[0])) / 2;
	//	const int32 h = int32(std::round(_ext[1])) / 2;
	//	const int32 cw = int32(x_ / 2);
	//	const int32 ch = int32(y_ / 2);

	//	//center kernel
	//	for (int32 x = cw - w; x <= cw + w; ++x) {
	//		for (int32 y = ch - h; y <= ch + h; ++y) {	
	//			const int32 i = kidx(x, y);
	//			fi_temp[i] = 1.f;
	//		}
	//	}
	//	/*
	//	for (int32 y = 0; y <= std::round(_ext[1]); ++y) {
	//		for (int32 x = 0; x <= std::round(_ext[0]); ++x) {
	//			const int32 i = idx(x, y);
	//			if (i < 0 || i >= map.size()) continue;
	//			fi_temp[i] = 1.f;
	//		}
	//	}
	//	*/

	//	if constexpr (std::is_same_v<float, T>)
	//		fftwf_execute_dft_r2c(std::get<0>(kplan_r2c), fi_temp.data(), reinterpret_cast<fftwf_complex*>(ci_temp.data()));
	//	else fftw_execute_dft_r2c(std::get<1>(kplan_r2c), fi_temp.data(), reinterpret_cast<fftw_complex*>(ci_temp.data()));
	//	
	//	for (int32 i = 0; i < ci_temp.size(); ++i)
	//		si_temp[i] = ci_temp[i];

	//	for (int32 x = 0; x < x_; ++x) {
	//		for (int32 y = 0; y < y_; ++y) {
	//			const int32 i1 = kidx(x, y);
	//			const int32 i2 = idx(x, y);
	//			si_temp[i2] = c_temp[i2] * ci_temp[i1];
	//		}
	//	}

	//	//for (int32 x = 0; x < x_; ++x) {
	//	//	for (int32 y = 0; y < y_; ++y) {
	//	//		const int32 i1 = kidx(x, y);
	//	//		result[1][i1] = int32(ci_temp[i1].real());
	//	//	}
	//	//}

	//	//return result;

	//	if constexpr (std::is_same_v<float, T>)
	//		fftwf_execute_dft_c2r(std::get<0>(plan_c2r), reinterpret_cast<fftwf_complex*>(si_temp.data()), fi_temp.data());
	//	else fftw_execute_dft_c2r(std::get<1>(plan_c2r), reinterpret_cast<fftw_complex*>(si_temp.data()), fi_temp.data());

	//	std::fill(c_temp.begin(), c_temp.end(), 0);

	//	for (int32 x = 0; x < x_; ++x) {
	//		for (int32 y = 0; y < y_; ++y) {
	//			const int32 i1 = kidx(x, y);
	//			const int32 i2 = idx(x, y);
	//			if (i2 < 0 || i2 >= c_temp.size()) continue;
	//			const int32 expec = int32(std::round(_ext[0]) + 1) * int32(std::round(_ext[0]) + 1) * map[i1];
	//			const int32 r = int32(std::abs(std::round(fi_temp[i2] * SCALE)));
	//			//result[0][i1] = r > expec ? -1 : r < expec ? 1 : 0;
	//			//result[1][i1] = r;
	//			c_temp[i2] = { std::abs(std::round(fi_temp[i2] * SCALE)) , 0. };
	//		}
	//	}

	//}

	//return c_temp;

	///*if (_perm | 0x2) {
	//}
	//
	//if (_perm | 0x8) {
	//}

	//if (_perm | 0x10) {
	//}

	//if (_perm | 0x20) {
	//}

	//if (_perm | 0x40) {
	//}*/

	//return result;
	
}//Util::HeightMap::overlap

template<class T>
void Util::HeightMap<T>::push(const FVector& _pos, const FVector& _ext) {

	const int32 minX = int32(_pos[0] - _ext[0]);
	const int32 maxX = int32(_pos[0] + _ext[0]);

	const int32 minY = int32(_pos[1] - _ext[1]);
	const int32 maxY = int32(_pos[1] + _ext[1]);

	const int32 maxZ = int32(_pos[2] + _ext[2]);

	for (int32 n0 = minY; n0 < maxY; ++n0) {
		for (int32 n1 = minX; n1 < maxX; ++n1) {
			const int32 i = n1 + n0 * n1_;
			map[i] = maxZ;
		}
	}

}//Util::HeightMap::push

template<class T>
void Util::HeightMap<T>::print_size_in_bytes() const {
	//std::cout << "Total size: " << double(map.size() * sizeof(int32) + f_temp.size() * sizeof(T) + c_temp.size() * sizeof(std::complex<T>) +
	//	fi_temp.size() * sizeof(T) + ci_temp.size() * sizeof(T) + 6 * result[0].size() * sizeof(int32)) / 1000000. << "Mb" << std::endl;
}//Util::HeightMap::print_size_in_bytes

template<class T>
int32 Util::HeightMap<T>::idx(const int32 _n0, const int32 _n1) const {
	const int32 o1 = p_n1_ / 2;
	const int32 o0 = p_n0_ / 2;

	const int32 out = (_n1 + o1) + (_n0 + o0) * p_n1_;
	/*if (out < 0 || out >= f_temp.size()) {
		std::cout << "idx: " << out << std::endl;
		return 0;
	}*/
	return out;
}//Util::HeightMap::idx

template<class T>
int32 Util::HeightMap<T>::kidx(const int32 _n0, const int32 _n1) const {
	const int32 out = _n1 + _n0 * n1_;
	if (out < 0 || out >= map.size()) {
		std::cout << "kidx: " << out << std::endl;
		return 0;
	}
	return out;
}//Util::HeightMap::idx
