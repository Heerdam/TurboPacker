#pragma once

#include "Util.h"

#include "GameFramework/Actor.h"
#include "lode/lodepng.h"

#include "DiscreteMap.generated.h"

//1200*800*1700 - 1,632,000,000
//  - 623,208,300

namespace TurboPacker {

	/*
		Z_: Up axis
		_n: negativ up axis
		XY_: object axis on the world axis
		_N: rotation of N*90° in ccw around the up axis
	*/
	UENUM()
	enum class EAxisPerm : uint8 {

		Z_XY_0, Z_XY_1, Z_XY_2, Z_XY_3,
		Y_XZ_0, Y_XZ_1, Y_XZ_2, Y_XZ_3,
		X_YZ_0, X_YZ_1, X_YZ_2, X_YZ_3,

		Z_n_XY_0, Z_n_XY_1, Z_n_XY_2, Z_n_XY_3,
		Y_n_XZ_0, Y_n_XZ_1, Y_n_XZ_2, Y_n_XZ_3,
		X_n_YZ_0, X_n_YZ_1, X_n_YZ_2, X_n_YZ_3

	};//EAxisPerm

	namespace Detail {
		FTransform make_transform(
			EAxisPerm _perm,
			const FBox& _target,
			const FVector& _pivot_offset,
			const FVector& _relative_offset
		);
	}

	namespace Spectral {

		namespace Detail {

			template<class T>
			struct FFTWDeleter;

			template<>
			struct FFTWDeleter<float> {
				void operator()(void* _ptr) const {
					fftwf_free(_ptr);
				}
			};

			template<>
			struct FFTWDeleter<std::complex<float>> {
				void operator()(void* _ptr) const {
					fftwf_free(_ptr);
				}
			};

			template<>
			struct FFTWDeleter<double> {
				void operator()(void* _ptr) const {
					fftw_free(_ptr);
				}
			};

			template<>
			struct FFTWDeleter<std::complex<double>> {
				void operator()(void* _ptr) const {
					fftw_free(_ptr);
				}
			};

			//-------------------

			template<class T>
			struct FFTPlanDeleter;

			template<>
			struct FFTPlanDeleter<float> {
				void operator()(void* _ptr) const {
					fftwf_destroy_plan((fftwf_plan)_ptr);
				}
			};

			template<>
			struct FFTPlanDeleter<double> {
				void operator()(void* _ptr) const {
					fftw_destroy_plan((fftw_plan)_ptr);
				}
			};

			//---------------------

			template <typename T>
			struct FFTWAlloc;

			template <>
			struct FFTWAlloc<float> {
				static void* malloc(int32 _size) {
					return fftwf_alloc_real(_size);
				}
			};

			template <>
			struct FFTWAlloc<double> {
				static void* malloc(int32 _size) {
					return fftw_alloc_real(_size);
				}
			};

			template <>
			struct FFTWAlloc<std::complex<float>> {
				static void* malloc(int32 _size) {
					return fftwf_alloc_complex(_size);
				}
			};

			template <>
			struct FFTWAlloc<std::complex<double>> {
				static void* malloc(int32 _size) {
					return fftw_alloc_complex(_size);
				}
			};

			//---------------------
			
			template <typename T>
			class FFTWAllocator {
			public:
				using value_type = T;

				FFTWAllocator() = default;

				template <class U>
				constexpr FFTWAllocator(const FFTWAllocator<U>&) noexcept {}

				T* allocate(std::size_t n) {
					return reinterpret_cast<T*>(FFTWAlloc<std::remove_pointer_t<T>>::malloc(n));
				}

				void deallocate(T* p, std::size_t n) noexcept {
					FFTWDeleter<std::remove_pointer_t<T>>()(p);
				}
			};

			template <class T, class U>
			bool operator==(const FFTWAllocator<T>&, const FFTWAllocator<U>&) { return true; }

			template <class T, class U>
			bool operator!=(const FFTWAllocator<T>&, const FFTWAllocator<U>&) { return false; }

			//---------------------

			template <typename T>
			struct PlanType;

			template <>
			struct PlanType<float> {
				using plan_t = std::remove_pointer_t<fftwf_plan>;
			};

			template <>
			struct PlanType<double> {
				using plan_t = std::remove_pointer_t<fftw_plan>;
			};

			template<class T>
			using PlantType_t = typename PlanType<T>::plan_t;

			//---------------------

			template<class T>
			struct FFTWPlaner;

			template <>
			struct FFTWPlaner<float> {
				static fftwf_plan r2c(int32 _n0, int32 _n1, float* _in, std::complex<float>* _out) {
					return fftwf_plan_dft_r2c_2d(_n0, _n1, _in, reinterpret_cast<fftwf_complex*>(_out), FFTW_MEASURE);
				}

				static fftwf_plan c2r(int32 _n0, int32 _n1, std::complex<float>* _in, float* _out) {
					return fftwf_plan_dft_c2r_2d(_n0, _n1, reinterpret_cast<fftwf_complex*>(_in), _out, FFTW_MEASURE);
				}
			};

			template <>
			struct FFTWPlaner<double> {
				static fftw_plan r2c(int32 _n0, int32 _n1, double* _in, std::complex<double>* _out) {
					return fftw_plan_dft_r2c_2d(_n0, _n1, _in, reinterpret_cast<fftw_complex*>(_out), FFTW_MEASURE);
				}

				static fftw_plan c2r(int32 _n0, int32 _n1, std::complex<double>* _in, double* _out) {
					return fftw_plan_dft_c2r_2d(_n0, _n1, reinterpret_cast<fftw_complex*>(_in), _out, FFTW_MEASURE);
				}
			};

			//---------------------

			template<class T>
			struct FFTWExecutor;

			template<>
			struct FFTWExecutor<float> {
				static void r2c(fftwf_plan _plan, float* _in, std::complex<float>* _out) {
					fftwf_execute_dft_r2c(_plan, _in, reinterpret_cast<fftwf_complex*>(_out));
				}

				static void c2r(fftwf_plan _plan, std::complex<float>* _in, float* _out) {
					fftwf_execute_dft_c2r(_plan, reinterpret_cast<fftwf_complex*>(_in), _out);
				}
			};

			template<>
			struct FFTWExecutor<double> {
				static void r2c(fftw_plan _plan, double* _in, std::complex<double>* _out) {
					fftw_execute_dft_r2c(_plan, _in, reinterpret_cast<fftw_complex*>(_out));
				}

				static void c2r(fftw_plan _plan, std::complex<double>* _in, double* _out) {
					fftw_execute_dft_c2r(_plan, reinterpret_cast<fftw_complex*>(_in), _out);
				}
			};

		}//Detail

		template<class T>
		using FFTWVector = std::vector<T, Detail::FFTWAllocator<T>>;

		template<class T>
		using FFTWUniquePtr = std::unique_ptr<T, Detail::FFTWDeleter<T>>;

		template<class T>
		using FFTWPlanPtr = std::unique_ptr<typename Detail::PlantType_t<T>, Detail::FFTPlanDeleter<T>>;

		namespace Debug {

			template<class C>
			void image_complex(
				const int32 _n0,
				const int32 _n1,
				C _buffer,
				const std::string& _filename,
				double _max = 2.0
			);

			template<class R, bool LogScaling>
			void image_real(
				const int32 _n0,
				const int32 _n1,
				R* _buffer,
				const std::string& _filename
			);

		}//Debug

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
			FFTWVector<T> signal_r;
			FFTWVector<std::complex<T>> signal_c;

			//kernel
			FFTWVector<T> kernel_r;
			FFTWVector<std::complex<T>> kernel_c;

			//sobel
			FFTWVector<std::complex<T>> sobel_x_c;
			FFTWVector<std::complex<T>> sobel_y_c;

			//gradient kernel
			FFTWVector<T> gradient_kernel_r ;
			FFTWVector<std::complex<T>> gradient_kernel_c;

			//gradient
			FFTWVector<T> gradient_r;
			FFTWVector<std::complex<T>> gradient_c;

			//temp
			FFTWVector<std::complex<T>> temp_c;

			//result
			mutable std::array<std::vector<char>, 6> result;

			//plans
			FFTWPlanPtr<T> plan_r2c;
			FFTWPlanPtr<T> plan_c2r;

		public:
			HeightMap(const int32 _n0, const int32 _n1, const int32 _h);
			HeightMap(HeightMap&&) = default;
			HeightMap(const HeightMap&) = delete;
			//--------------------
			// Z_XY | ZYX | Y_XZ | Y_ZX | X_YZ | X_ZY
			template<bool DEBUG = false>
			const std::array<std::vector<char>, 6>& overlap(const FVector& _ext, const int32 _perm = 0x3F);
			void push(const FVector& _pos, const FVector& _ext);
			void print_size_in_bytes() const;

		};//HeightMap

	}//Spectral

}//TurboPacker

UCLASS(Blueprintable, BlueprintType)
class TURBOPACKER_API ASpectralPacker : public AActor {
	GENERATED_BODY()

public:

	UPROPERTY(EditAnywhere)
	FIntVector Bounds;

	UPROPERTY(EditDefaultsOnly)
	TArray<TSubclassOf<class APackerBox>> Boxes;

	UFUNCTION(BlueprintCallable, CallInEditor, Category = Packer)
	void Clear();

	UFUNCTION(BlueprintCallable, CallInEditor, Category = Packer)
	void Pack();

};//ASpectralPacker

//------------------------------------------

UCLASS(Blueprintable, BlueprintType)
class TURBOPACKER_API ASpectralTester : public AActor {
	GENERATED_BODY()

public:
	UFUNCTION(BlueprintCallable, CallInEditor, Category = Test)
	void Clear();

	UFUNCTION(BlueprintCallable, CallInEditor, Category = Test)
	void TestHeightMap();

	UFUNCTION(BlueprintCallable, CallInEditor, Category = Test)
	void TestFFTW();

	UFUNCTION(BlueprintCallable, CallInEditor, Category = Test)
	void TestKernel();

	UPROPERTY(EditAnywhere)
	int32 Itensity = 8;

};//ASpectralTester

UCLASS(Blueprintable, BlueprintType)
class TURBOPACKER_API APackerBox : public AActor {
	GENERATED_BODY()

public:

};//APackerBox

//------------------------------------------------
//------------------ HeightMap -------------------
//------------------------------------------------

template<class T>
TurboPacker::Spectral::HeightMap<T>::HeightMap(const int32 _n0, const int32 _n1, const int32 _h) :
	n0_(_n0), n1_(_n1), n1_c_(_n1 / 2 + 1), h_(_h),
	p_n0_(2 * _n0 - 1), p_n1_(2 * _n1 - 1), p_n1_c_(p_n1_ / 2 + 1) {

	using namespace Util;
	using namespace Detail;

	map.resize(_n0 * _n1);
	std::fill(map.begin(), map.end(), 0);
	for (int32 i = 0; i < 6; ++i)
		result[i].resize(_n0 * _n1);
	//-------------------------
	signal_r.resize(p_n0_ * p_n1_);
	signal_c.resize(p_n0_ * p_n1_);

	kernel_r.resize(p_n0_ * p_n1_);
	kernel_c.resize(p_n0_ * p_n1_);

	sobel_x_c.resize(p_n0_ * p_n1_);
	sobel_y_c.resize(p_n0_ * p_n1_);

	gradient_kernel_r.resize(p_n0_ * p_n1_);
	gradient_kernel_c.resize(p_n0_ * p_n1_);

	gradient_r.resize(p_n0_ * p_n1_);
	gradient_c.resize(p_n0_ * p_n1_);

	temp_c.resize(p_n0_ * p_n1_);

	plan_r2c = FFTWPlanPtr<T>(FFTWPlaner<T>::r2c(p_n0_, p_n1_, signal_r.data(), signal_c.data()));
	plan_c2r = FFTWPlanPtr<T>(FFTWPlaner<T>::c2r(p_n0_, p_n1_, signal_c.data(), signal_r.data()));
	//---------------
	FFTWVector<T> sobel_r(p_n0_ * p_n1_);
	sobel_r[0 + 0 * p_n1_] = -1.;
	sobel_r[1 + 0 * p_n1_] = 0.;
	sobel_r[2 + 0 * p_n1_] = 1.;
	sobel_r[0 + 1 * p_n1_] = -2.;
	sobel_r[1 + 1 * p_n1_] = 0.;
	sobel_r[2 + 1 * p_n1_] = 2.;
	sobel_r[0 + 2 * p_n1_] = -1.;
	sobel_r[1 + 2 * p_n1_] = 0.;
	sobel_r[2 + 2 * p_n1_] = 1.;
	FFTWExecutor<T>::r2c(plan_r2c.get(), sobel_r.data(), sobel_x_c.data());
	//---------------
	sobel_r[0 + 0 * p_n1_] = -1.;
	sobel_r[1 + 0 * p_n1_] = -2.;
	sobel_r[2 + 0 * p_n1_] = -1.;
	sobel_r[0 + 1 * p_n1_] = 0.;
	sobel_r[1 + 1 * p_n1_] = 0.;
	sobel_r[2 + 1 * p_n1_] = 0.;
	sobel_r[0 + 2 * p_n1_] = 1.;
	sobel_r[1 + 2 * p_n1_] = 2.;
	sobel_r[2 + 2 * p_n1_] = 1.;
	FFTWExecutor<T>::r2c(plan_r2c.get(), sobel_r.data(), sobel_y_c.data());
}//TurboPacker::Spectral::HeightMap::HeightMap

template<class T>
template<bool DEBUG>
const std::array<std::vector<char>, 6>& TurboPacker::Spectral::HeightMap<T>::overlap(const FVector& _ext, const int32 _perm) {
	using namespace Util;
	using namespace Debug;
	using namespace Detail;
	//----------
	const int32 o1 = n1_ / 2;
	const int32 o0 = n0_ / 2;
	//------------ Signal -----------
	for (int32 i = 0; i < p_n0_ * p_n1_; ++i)
		signal_r[i] = 0.;
	for (int32 n0 = 0; n0 < n0_; ++n0) {
		for (int32 n1 = 0; n1 < n1_; ++n1) {	
			const int32 i1 = n1 + n0 * n1_;
			const int32 i2 = (n1 + o1) + (n0 + o0) * p_n1_;
			signal_r[i2] = T(map[i1]);
		}
	}

	if constexpr (DEBUG) image_real<T, false>(p_n0_, p_n1_, signal_r.data(), "1_signal.png");

	FFTWExecutor<T>::r2c(plan_r2c.get(), signal_r.data(), signal_c.data());

	if constexpr (DEBUG) image_complex(p_n0_, p_n1_c_, signal_c.data(), "2_signal.png", 2.0);

	if (_perm | 0x1) {
		const int32 w = int32(std::floor(_ext[0]));
		const int32 h = int32(std::floor(_ext[1]));

		//------------ Kernel -----------
		//needs to be in the corner
		for (int32 i = 0; i < p_n0_ * p_n1_; ++i)
			kernel_r[i] = 0.;

		for (int32 n0 = 0; n0 < h; ++n0) {
			for (int32 n1 = 0; n1 < w; ++n1) {
				const int32 i = n1 + n0 * p_n1_;
				kernel_r[i] = 1.f;
			}
		}

		if constexpr (DEBUG) image_real<T, false>(p_n0_, p_n1_, kernel_r.data(), "3_kernel.png");

		FFTWExecutor<T>::r2c(plan_r2c.get(), kernel_r.data(), kernel_c.data());

		if constexpr (DEBUG) image_complex(p_n0_, p_n1_c_, kernel_c.data(), "4_kernel.png", 2.0);

		//------------ Convolution -----------
		for (int32 n0 = 0; n0 < p_n0_; ++n0) {
			for (int32 n1 = 0; n1 < p_n1_c_; ++n1) {
				const int32 i = n1 + n0 * p_n1_c_;
				const std::complex<T>& c1 = signal_c[i];
				const std::complex<T>& c2 = kernel_c[i];
				const auto r = c1 * c2;
				temp_c[i] = r;
			}
		}

		FFTWExecutor<T>::c2r(plan_c2r.get(), temp_c.data(), signal_r.data());

		if constexpr (DEBUG) image_real<T, false>(p_n0_, p_n1_, signal_r.data(), "5_convolution.png");

		const int32 w2 = w / 2;
		const int32 h2 = h / 2;
		const T SCALE = 1. / T(p_n0_ * p_n1_);

		for (int32 n0 = 0; n0 < n0_; ++n0) {
			for (int32 n1 = 0; n1 < n1_; ++n1) {
				const int32 i1 = n1 + n0 * n1_;
				const int32 i2 = (n1 + o1 + h2) + (n0 + o0 + w2) * p_n1_;

				const int32 expec = w * h * map[i1];
				const int32 r = int32(std::abs(std::round(signal_r[i2] * SCALE)));
				signal_r[i2] = r > expec ? -1. : r < expec ? 1. : 0.;
				//result[_idx][i1] = r > expec ? -1 : r < expec ? 1 : 0;
			}
		}

		if constexpr (DEBUG) image_real<T, false>(p_n0_, p_n1_, signal_r.data(), "6_convolution.png");

		/*
		//border kernel
		for (int32 n0 = 0; n0 < p_n0_; ++n0) {
			for (int32 n1 = 0; n1 < p_n1_; ++n1) {
				const int32 i = n1 + n0 * p_n1_;

				const int32 v = int32(std::round(kernel_r[i]));
				if (v == 0) continue;

				const int32 i1 = (n1 + 1) + n0 * p_n1_;
				const int32 i2 = n1 + (n0 - 1) * p_n1_;
				const int32 i3 = (n1 - 1) + n0 * p_n1_;
				const int32 i4 = n1 + (n0 + 1) * p_n1_;

				gradient_kernel_r[i1] = 1.;
				gradient_kernel_r[i2] = 1.;
				gradient_kernel_r[i3] = 1.;
				gradient_kernel_r[i4] = 1.;
			}
		}

		for (int32 n0 = 0; n0 < p_n0_; ++n0) {
			for (int32 n1 = 0; n1 < p_n1_; ++n1) {
				const int32 i = n1 + n0 * p_n1_;

				const int32 v = int32(std::round(kernel_r[i]));
				if (v == 0) continue;


			}
		}
		*/

		//gradient
		for (int32 n0 = 0; n0 < p_n0_; ++n0) {
			for (int32 n1 = 0; n1 < p_n1_; ++n1) {
				const int32 i = n1 + n0 * p_n1_;

				const int32 v = int32(std::round(signal_r[i]));
				if (v == 0) continue;

				const int32 i1 = (n1 + 1) + n0 * p_n1_;
				const int32 i2 = n1 + (n0 - 1) * p_n1_;
				const int32 i3 = (n1 - 1) + n0 * p_n1_;
				const int32 i4 = n1 + (n0 + 1) * p_n1_;

				const T v1 = int32(std::round(signal_r[i1]));
				const T v2 = int32(std::round(signal_r[i2]));
				const T v3 = int32(std::round(signal_r[i3]));
				const T v4 = int32(std::round(signal_r[i4]));

				if (v1 != 0 && v != v1) {
					gradient_r[i] = 1.;
					continue;
				}
				
				if (v2 != 0 && v != v2) {
					gradient_r[i] = 1.;
					continue;
				}

				if (v3 != 0 && v != v3) {
					gradient_r[i] = 1.;
					continue;
				}

				if (v4 != 0 && v != v4) {
					gradient_r[i] = 1.;
					continue;
				}

				gradient_r[i] = 0.;

			}
		}

		if constexpr (DEBUG) image_real<T, false>(p_n0_, p_n1_, gradient_r.data(), "8_gradient.png");

		//convolution

		//result

	}

	/*
	if (_perm | 0x2) {
		const int32 w = int32(std::floor(_ext[0]));
		const int32 h = int32(std::floor(_ext[1]));
		overlap_int(1, h, w);
	}

	if (_perm | 0x8) {
		const int32 w = int32(std::floor(_ext[0]));
		const int32 h = int32(std::floor(_ext[2]));
		overlap_int(2, w, h);
	}

	if (_perm | 0x10) {
		const int32 w = int32(std::floor(_ext[0]));
		const int32 h = int32(std::floor(_ext[2]));
		overlap_int(3, h, w);
	}

	if (_perm | 0x20) {
		const int32 w = int32(std::floor(_ext[1]));
		const int32 h = int32(std::floor(_ext[2]));
		overlap_int(4, w, h);
	}

	if (_perm | 0x40) {
		const int32 w = int32(std::floor(_ext[1]));
		const int32 h = int32(std::floor(_ext[2]));
		overlap_int(5, h, w);
	}
	*/

	return result;

}//TurboPacker::Spectral::HeightMap::overlap

template<class T>
void TurboPacker::Spectral::HeightMap<T>::push(const FVector& _pos, const FVector& _ext) {

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

}//TurboPacker::Spectral::HeightMap::push

template<class T>
void TurboPacker::Spectral::HeightMap<T>::print_size_in_bytes() const {
	//std::cout << "Total size: " << double(map.size() * sizeof(int32) + f_temp.size() * sizeof(T) + c_temp.size() * sizeof(std::complex<T>) +
	//	fi_temp.size() * sizeof(T) + ci_temp.size() * sizeof(T) + 6 * result[0].size() * sizeof(int32)) / 1000000. << "Mb" << std::endl;
}//TurboPacker::Spectral::HeightMap::print_size_in_bytes

//----------------------------------------

template<class C>
void TurboPacker::Spectral::Debug::image_complex(
	const int32 _n0,
	const int32 _n1,
	C _buffer,
	const std::string& _filename,
	double _max
) {
	{
		const std::filesystem::path path = std::filesystem::path(TCHAR_TO_UTF8(*FPaths::ProjectDir())) / "overlap";
		if (!std::filesystem::exists(path))
			std::filesystem::create_directory(path);
	}

	std::vector<unsigned char> img;
	for (int32 n0 = 0; n0 < _n0; ++n0) {
		for (int32 n1 = 0; n1 < _n1; ++n1) {
			const int32 i1 = n1 + n0 * _n1;

			const std::complex<double> c(_buffer[i1].real(), _buffer[i1].imag());
			const FVector col = Util::HSL_to_RGB_rad<double>(Util::c_to_HSL<double>(c, _max));

			img.push_back(255 - (unsigned char)(col[0] * 255));
			img.push_back(255 - (unsigned char)(col[1] * 255));
			img.push_back(255 - (unsigned char)(col[2] * 255));
			img.push_back(255);
		}
	}
	const std::filesystem::path path = std::filesystem::path(TCHAR_TO_UTF8(*FPaths::ProjectDir())) / "overlap" / _filename;
	lodepng::encode(path.string(), img.data(), _n1, _n0);

}//TurboPacker::Spectral::Debug::image_complex

template<class R, bool LogScaling>
void TurboPacker::Spectral::Debug::image_real(
	const int32 _n0,
	const int32 _n1,
	R* _buffer,
	const std::string& _filename
) {
	{
		const std::filesystem::path path = std::filesystem::path(TCHAR_TO_UTF8(*FPaths::ProjectDir())) / "overlap";
		if (!std::filesystem::exists(path))
			std::filesystem::create_directory(path);
	}

	double max = -std::numeric_limits<double>::infinity();
	for (int32 i = 0; i < _n0 * _n1; ++i)
		max = std::max<double>(max, std::abs(_buffer[i]));

	std::vector<unsigned char> img;
	for (int32 n0 = 0; n0 < _n0; ++n0) {
		for (int32 n1 = 0; n1 < _n1; ++n1) {
			const int32 i1 = n1 + n0 * _n1;
			unsigned char c = 0;

			if constexpr (LogScaling) {
				c = std::abs(1. - std::abs(std::log(std::abs(1. + _buffer[i1]))) / std::log(max));
			} else {
				const double frac = 1. / max;
				c = 255 * (_buffer[i1] * frac);
			}

			img.push_back(c);
			img.push_back(c);
			img.push_back(c);
			img.push_back(255);
		}
	}
	const std::filesystem::path path = std::filesystem::path(TCHAR_TO_UTF8(*FPaths::ProjectDir())) / "overlap/" / _filename;
	lodepng::encode(path.string(), img.data(), _n1, _n0);
}
