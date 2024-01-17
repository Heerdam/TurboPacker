#pragma once

#include "Util.h"

#include "GameFramework/Actor.h"
//#include "lode/lodepng.h"
//#include <immintrin.h>
//
//#include "DiscreteMap.generated.h"
//
////1200*800*1700 - 1,632,000,000
////  - 623,208,300
//
//namespace TurboPacker {
//
//	/*
//		Z_: Up axis
//		_n: negativ up axis
//		XY_: object axis on the world axis
//		_N: rotation of N*90° in ccw around the up axis
//	*/
//	UENUM()
//	enum class EAxisPerm : uint8 {
//
//		Z_XY_0, Z_XY_1, Z_XY_2, Z_XY_3,
//		Y_XZ_0, Y_XZ_1, Y_XZ_2, Y_XZ_3,
//		X_YZ_0, X_YZ_1, X_YZ_2, X_YZ_3,
//
//		Z_n_XY_0, Z_n_XY_1, Z_n_XY_2, Z_n_XY_3,
//		Y_n_XZ_0, Y_n_XZ_1, Y_n_XZ_2, Y_n_XZ_3,
//		X_n_YZ_0, X_n_YZ_1, X_n_YZ_2, X_n_YZ_3
//
//	};//EAxisPerm
//
//	namespace Detail {
//		FTransform make_transform(
//			EAxisPerm _perm,
//			const FBox& _target,
//			const FVector& _pivot_offset,
//			const FVector& _relative_offset
//		);
//	}
//
//	namespace Spectral {
//
//		namespace Detail {
//
//			template<class T>
//			struct FFTWDeleter;
//
//			template<>
//			struct FFTWDeleter<float> {
//				void operator()(void* _ptr) const {
//					fftwf_free(_ptr);
//				}
//			};
//
//			template<>
//			struct FFTWDeleter<std::complex<float>> {
//				void operator()(void* _ptr) const {
//					fftwf_free(_ptr);
//				}
//			};
//
//			template<>
//			struct FFTWDeleter<double> {
//				void operator()(void* _ptr) const {
//					fftw_free(_ptr);
//				}
//			};
//
//			template<>
//			struct FFTWDeleter<std::complex<double>> {
//				void operator()(void* _ptr) const {
//					fftw_free(_ptr);
//				}
//			};
//
//			//-------------------
//
//			template<class T>
//			struct FFTPlanDeleter;
//
//			template<>
//			struct FFTPlanDeleter<float> {
//				void operator()(void* _ptr) const {
//					fftwf_destroy_plan((fftwf_plan)_ptr);
//				}
//			};
//
//			template<>
//			struct FFTPlanDeleter<double> {
//				void operator()(void* _ptr) const {
//					fftw_destroy_plan((fftw_plan)_ptr);
//				}
//			};
//
//			//---------------------
//
//			template <typename T>
//			struct FFTWAlloc;
//
//			template <>
//			struct FFTWAlloc<float> {
//				static void* malloc(int32 _size) {
//					return fftwf_alloc_real(_size);
//				}
//			};
//
//			template <>
//			struct FFTWAlloc<double> {
//				static void* malloc(int32 _size) {
//					return fftw_alloc_real(_size);
//				}
//			};
//
//			template <>
//			struct FFTWAlloc<std::complex<float>> {
//				static void* malloc(int32 _size) {
//					return fftwf_alloc_complex(_size);
//				}
//			};
//
//			template <>
//			struct FFTWAlloc<std::complex<double>> {
//				static void* malloc(int32 _size) {
//					return fftw_alloc_complex(_size);
//				}
//			};
//
//			//---------------------
//			
//			template <class T>
//			class FFTWAllocator {
//			public:
//				using value_type = T;
//
//				FFTWAllocator() = default;
//
//				template <class U>
//				constexpr FFTWAllocator(const FFTWAllocator<U>&) noexcept {}
//
//				T* allocate(std::size_t n) {
//					return reinterpret_cast<T*>(FFTWAlloc<std::remove_pointer_t<T>>::malloc(n));
//				}
//
//				void deallocate(T* p, std::size_t n) noexcept {
//					FFTWDeleter<std::remove_pointer_t<T>>()(p);
//				}
//			};
//
//			template <class T, class U>
//			bool operator==(const FFTWAllocator<T>&, const FFTWAllocator<U>&) { return true; }
//
//			template <class T, class U>
//			bool operator!=(const FFTWAllocator<T>&, const FFTWAllocator<U>&) { return false; }
//
//			//---------------------
//
//			template <typename T>
//			struct PlanType;
//
//			template <>
//			struct PlanType<float> {
//				using plan_t = std::remove_pointer_t<fftwf_plan>;
//			};
//
//			template <>
//			struct PlanType<double> {
//				using plan_t = std::remove_pointer_t<fftw_plan>;
//			};
//
//			template<class T>
//			using PlantType_t = typename PlanType<T>::plan_t;
//
//			//---------------------
//
//			template<class T>
//			struct FFTWPlaner;
//
//			template <>
//			struct FFTWPlaner<float> {
//				static fftwf_plan r2c(int32 _n0, int32 _n1, float* _in, std::complex<float>* _out) {
//					return fftwf_plan_dft_r2c_2d(_n0, _n1, _in, reinterpret_cast<fftwf_complex*>(_out), FFTW_MEASURE);
//				}
//
//				static fftwf_plan c2r(int32 _n0, int32 _n1, std::complex<float>* _in, float* _out) {
//					return fftwf_plan_dft_c2r_2d(_n0, _n1, reinterpret_cast<fftwf_complex*>(_in), _out, FFTW_MEASURE);
//				}
//			};
//
//			template <>
//			struct FFTWPlaner<double> {
//				static fftw_plan r2c(int32 _n0, int32 _n1, double* _in, std::complex<double>* _out) {
//					return fftw_plan_dft_r2c_2d(_n0, _n1, _in, reinterpret_cast<fftw_complex*>(_out), FFTW_MEASURE);
//				}
//
//				static fftw_plan c2r(int32 _n0, int32 _n1, std::complex<double>* _in, double* _out) {
//					return fftw_plan_dft_c2r_2d(_n0, _n1, reinterpret_cast<fftw_complex*>(_in), _out, FFTW_MEASURE);
//				}
//			};
//
//			//---------------------
//
//			template<class T>
//			struct FFTWExecutor;
//
//			template<>
//			struct FFTWExecutor<float> {
//				static void r2c(fftwf_plan _plan, float* _in, std::complex<float>* _out) {
//					fftwf_execute_dft_r2c(_plan, _in, reinterpret_cast<fftwf_complex*>(_out));
//				}
//
//				static void c2r(fftwf_plan _plan, std::complex<float>* _in, float* _out) {
//					fftwf_execute_dft_c2r(_plan, reinterpret_cast<fftwf_complex*>(_in), _out);
//				}
//			};
//
//			template<>
//			struct FFTWExecutor<double> {
//				static void r2c(fftw_plan _plan, double* _in, std::complex<double>* _out) {
//					fftw_execute_dft_r2c(_plan, _in, reinterpret_cast<fftw_complex*>(_out));
//				}
//
//				static void c2r(fftw_plan _plan, std::complex<double>* _in, double* _out) {
//					fftw_execute_dft_c2r(_plan, reinterpret_cast<fftw_complex*>(_in), _out);
//				}
//			};
//
//			//---------------------
//
//			inline int32 power_of_2(int32 _n) {
//				if (_n == 0) return 1;
//				if ((_n & (_n - 1)) == 0) return _n;
//				int32 res = 1;
//				while (res < _n)
//					res *= 2;
//				return res;
//			};//power_of_2
//
//		}//Detail
//
//		template<class T>
//		using FFTWVector = std::vector<T, Detail::FFTWAllocator<T>>;
//
//		template<class T>
//		using FFTWUniquePtr = std::unique_ptr<T, Detail::FFTWDeleter<T>>;
//
//		template<class T>
//		using FFTWPlanPtr = std::unique_ptr<typename Detail::PlantType_t<T>, Detail::FFTPlanDeleter<T>>;
//
//		namespace Debug {
//
//			template<class C>
//			inline void image_complex(
//				const int32 _n0,
//				const int32 _n1,
//				C _buffer,
//				const std::string& _filename,
//				double _max = 2.0
//			);
//
//			template<class R, bool LogScaling>
//			inline void image_real(
//				const int32 _n0,
//				const int32 _n1,
//				R* _buffer,
//				const std::string& _filename
//			);
//
//		}//Debug
//
//		namespace Impl {
//
//			template<class T_, bool SIMD_, bool DEBUG_>
//			struct Config {
//				using T = T_;
//				constexpr static bool SIMD = SIMD_;
//				constexpr static bool DEBUG = DEBUG_;
//				//---------------
//				const int32 n0_;
//				const int32 n1_;
//				const int32 p_n0_;
//				const int32 p_n1_;
//				const int32 p_n1_c_;
//				const int32 h_;
//				//---------------
//				const int32 off0_;
//				const int32 off1_;
//				const T SCALE;
//				//---------------
//				UWorld* world = nullptr;
//				//---------------
//				FFTWVector<std::complex<T>> sobel_x_c;
//				FFTWVector<std::complex<T>> sobel_y_c;
//				//---------------
//				FFTWVector<T> temp1_r;
//				FFTWVector<T> temp2_r;
//				FFTWVector<std::complex<T>> temp_c;
//				//---------------
//				FFTWPlanPtr<T> plan_r2c;
//				FFTWPlanPtr<T> plan_c2r;
//				//---------------
//				Config() = delete;
//				Config(const Config&) = delete;
//				Config(Config&&) = default;
//				Config(int32 _n0, int32 _n1, int32 _h);
//				//---------------
//				std::ostream& operator<<(std::ostream& s);
//			};//Config
//
//			template<class T_, bool SIMD_, bool DEBUG_>
//			void fft_signal(
//				Config<T_, SIMD_, DEBUG_>& _config,
//				const FFTWVector<T_>& _map,
//				const FFTWVector<T_>& _map_log,
//				FFTWVector<std::complex<T_>>& _fft_signal_out,
//				FFTWVector<std::complex<T_>>& _fft_signal_log_out
//			);
//
//			template<class T_, class VOXLER_, bool SIMD_, bool DEBUG_>
//			void fft_kernel(
//				Config<T_, SIMD_, DEBUG_>& _config,
//				const VOXLER_& _voxler,
//				FFTWVector<std::complex<T_>>& _kernel_out,
//				FFTWVector<std::complex<T_>>& _border_kernel_out
//			);
//
//			template<class T_, bool SIMD_, bool DEBUG_>
//			void correlate_kernel_signal(
//				Config<T_, SIMD_, DEBUG_>& _config,
//				const FFTWVector<std::complex<T_>>& _kernel_in,
//				const FFTWVector<std::complex<T_>>& _signal_in,
//				const FFTWVector<std::complex<T_>>& _signal_log_in,
//				FFTWVector<T_>& _corr_out,
//				FFTWVector<T_>& _corr_log_out
//			);
//
//			template<class T_, bool SIMD_, bool DEBUG_>
//			void extract_solution(
//				Config<T_, SIMD_, DEBUG_>& _config,
//				const FFTWVector<T_>& _corr_in,
//				const FFTWVector<T_>& _corr_log_in,
//				FFTWVector<T_>& _sol_out
//			);
//
//			template<class T_, bool SIMD_, bool DEBUG_>
//			void gradients_signal_stencil(
//				Config<T_, SIMD_, DEBUG_>& _config,
//				const FFTWVector<T_>& _sol,
//				FFTWVector<std::complex<T_>>& _border_out
//			);
//
//			template<class T_, bool SIMD_, bool DEBUG_>
//			void gradients_signal(
//				Config<T_, SIMD_, DEBUG_>& _config,
//				const FFTWVector<std::complex<T_>>& _signal_in,
//				FFTWVector<T_>& _border_out
//			);
//
//			template<class T_, bool SIMD_, bool DEBUG_>
//			void proximity(
//				Config<T_, SIMD_, DEBUG_>& _config,
//				const FFTWVector<T_>& _border,
//				const FFTWVector<std::complex<T_>>& _border_kernel,
//				FFTWVector<T_>& _solution
//			);
//
//		}//Impl
//
//		//---------------------------------
//
//		template<class T, bool SIMD, bool DEBUG>
//		class TURBOPACKER_API MedianQuadTree {
//
//			struct Bucket {
//				//------------------
//				MedianQuadTree* p = nullptr;
//				FIntVector2 bmin, bmax;
//				//------------------
//				const int32 level;
//				T median, max, min;
//				bool isFlat = false;
//				std::vector<FIntVector2> l;
//				std::vector<FIntVector2> h;
//				Bucket(
//					MedianQuadTree* _p, 
//					const FIntVector2& _bmin,
//					const FIntVector2& _bmax,
//					int32 _l
//				) : 
//					p(_p), bmin(_bmin), bmax(_bmax), level(_l) {}
//				void recompute();
//				std::pair<int32, int32> overlap(const FIntVector2& _min, const FIntVector2& _max, T _h) const;
//				//------------------
//				std::ostream& operator<< (std::ostream& s) const;
//			};//Bucket
//
//			struct Node {
//				//------------------
//				MedianQuadTree* p = nullptr;
//				FIntVector2 bmin, bmax;
//				//------------------
//				const int32 level;
//				T median, max, min;
//				bool isFlat = false;
//				std::array<int32, 2> l, h;
//				std::array<std::variant<std::unique_ptr<Node>, std::unique_ptr<Bucket>>, 4> c;
//				Node(
//					MedianQuadTree* _p,
//					const FIntVector2& _bmin,
//					const FIntVector2& _bmax,
//					int32 _l
//				) :
//					p(_p), bmin(_bmin), bmax(_bmax), level(_l) {}
//				void recompute();
//				std::pair<int32, int32> overlap(const FIntVector2& _min, const FIntVector2& _max, T _h) const;
//				//------------------
//				std::ostream& operator<< (std::ostream& s) const;
//			};//Node
//
//			//friend Bucket;
//
//			std::unique_ptr<Node> root;
//
//			const FFTWVector<T>& map;
//			const int32 n0_;
//			const int32 n1_;
//
//			const int32 extend;
//			int32 maxLevel;
//
//			UWorld* world = nullptr;
//
//		public:
//			using TYPE = T;
//			constexpr static bool USE_SIMD = SIMD;
//			constexpr static bool USE_DEBUG = DEBUG;
//			//----------------
//			MedianQuadTree(
//				Impl::Config<T, SIMD, DEBUG>& _config,
//				const FFTWVector<T>& _map,
//				const int32 _min_ext = 12
//			);
//
//			void recompute();
//			std::pair<int32, int32> check_height(const FIntVector& _pos, const FIntVector2& _ext) const;
//			std::ostream& operator<< (std::ostream& s) const;
//		};//MedianQuadTree
//
//		//---------------------------------
//
//		template<class T>
//		class TURBOPACKER_API BoxStack {
//
//			int32 idd = 0;
//			std::vector<std::tuple<int32, int32, int32>> list;
//			tsl::robin_map<int32, UE::Math::TBox<T>> boxes;
//
//			std::vector<UE::Math::TBox<T>> boxs;
//
//		public:
//
//			void push(const UE::Math::TBox<T>& _box);
//			bool overlap(const UE::Math::TBox<T>& _box) const;
//			bool overlap_naive(const UE::Math::TBox<T>& _box) const;
//
//		};//BoxStack
//
//		//---------------------------------
//
//		template<class T, bool SIMD, bool DEBUG>
//		class TURBOPACKER_API HeightMap {
//
//			static_assert(sizeof(std::complex<float>) == sizeof(fftwf_complex));
//			static_assert(sizeof(std::complex<double>) == sizeof(fftw_complex));
//			static_assert(std::is_same_v<T, double> || std::is_same_v<T, float>);
//
//		public:
//			using TYPE = T;
//			constexpr static bool USE_SIMD = SIMD;
//			constexpr static bool USE_DEBUG = DEBUG;
//
//		private:
//
//			std::unique_ptr<Impl::Config<T, SIMD, DEBUG>> config_;
//
//			//map
//			FFTWVector<T> map;
//			FFTWVector<T> map_log;
//
//			//signal
//			FFTWVector<std::complex<T>> signal_c;
//			FFTWVector<std::complex<T>> signal_log_c;
//
//			//kernel
//			FFTWVector<std::complex<T>> kernel_c;
//
//			//sobel
//			FFTWVector<std::complex<T>> border_c;
//
//			//plans
//			FFTWPlanPtr<T> plan_r2c;
//			FFTWPlanPtr<T> plan_c2r;
//
//			//------------------------
//
//
//
//			//n0, n1, h, c
//			std::function<T(int32, int32, int32, int32)> costf;
//
//		public:
//			HeightMap(const int32 _n0, const int32 _n1, const int32 _h);
//			HeightMap(HeightMap&&) = default;
//			HeightMap(const HeightMap&) = delete;
//			//--------------------
//			// Z_XY | ZYX | Y_XZ | Y_ZX | X_YZ | X_ZY
//			std::vector<std::tuple<FIntVector, T, int32>> overlap(const FVector& _ext, const int32 _perm = 0x3F);
//			void push(const FVector& _pos, const FVector& _ext);
//			void print_size_in_bytes() const;
//			void save_heightmap(int32 _id);
//			void set_cost_function(const std::function<T(int32, int32, int32, int32)>& _cf) { costf = _cf; }
//		};//HeightMap
//
//	}//Spectral
//
//}//TurboPacker
//
//UCLASS(Blueprintable, BlueprintType)
//class TURBOPACKER_API ASpectralPacker : public AActor {
//	GENERATED_BODY()
//
//public:
//
//	UPROPERTY(EditAnywhere)
//	FIntVector Bounds;
//
//	UPROPERTY(EditDefaultsOnly)
//	TArray<TSubclassOf<class APackerBox>> Boxes;
//
//	UFUNCTION(BlueprintCallable, CallInEditor, Category = Packer)
//	void Clear();
//
//	UFUNCTION(BlueprintCallable, CallInEditor, Category = Packer)
//	void Pack();
//
//};//ASpectralPacker
//
////------------------------------------------
//
//UCLASS(Blueprintable, BlueprintType)
//class TURBOPACKER_API ASpectralTester : public AActor {
//	GENERATED_BODY()
//
//public:
//
//	UPROPERTY(EditDefaultsOnly)
//	TArray<TSubclassOf<class APackerBox>> Boxes;
//
//	UPROPERTY(EditAnywhere)
//	FIntVector Bounds = FIntVector(1200, 800, 1700);
//
//	UFUNCTION(BlueprintCallable, CallInEditor, Category = Test)
//	void Clear();
//
//	UFUNCTION(BlueprintCallable, CallInEditor, Category = Test)
//	void TestPacker();
//
//	UFUNCTION(BlueprintCallable, CallInEditor, Category = Test)
//	void TestHeightMap();
//
//	UFUNCTION(BlueprintCallable, CallInEditor, Category = Test)
//	void TestKernel();
//
//	UFUNCTION(BlueprintCallable, CallInEditor, Category = Test)
//	void BenchConv();
//
//	UFUNCTION(BlueprintCallable, CallInEditor, Category = Test)
//	void QuadTreeTester();
//
//	UFUNCTION(BlueprintCallable, CallInEditor, Category = Test)
//	void QuadTreeTester2();
//
//	UFUNCTION(BlueprintCallable, CallInEditor, Category = Test)
//	void QuadTreeBench();
//
//	UFUNCTION(BlueprintCallable, CallInEditor, Category = Test)
//	void BoxStackBench();
//
//};//ASpectralTester
//
//UCLASS(Blueprintable, BlueprintType)
//class TURBOPACKER_API APackerBox : public AActor {
//	GENERATED_BODY()
//
//public:
//
//	UPROPERTY(EditDefaultsOnly, Category = Boxler)
//	UStaticMeshComponent* mesh = nullptr;
//
//	APackerBox();
//	FBox get_aabb() const;
//	virtual FVector get_relative_location() const { return FVector(0.); }
//
//};//APackerBox
//
////------------------------------------------------
////------------------- Config ---------------------
////------------------------------------------------
//
//template<class T_, bool SIMD_, bool DEBUG_>
//TurboPacker::Spectral::Impl::Config<T_, SIMD_, DEBUG_>::Config(
//	int32 _n0, 
//	int32 _n1,
//	int32 _h
//) :
//	n0_(_n0), n1_(_n1),
//	p_n0_(Detail::power_of_2(int32(std::pow(std::floor(std::sqrt(T(2 * _n0 - 1))) + 1, 2)))),
//	p_n1_(Detail::power_of_2(int32(std::pow(std::floor(std::sqrt(T(2 * _n1 - 1))) + 1, 2)))),
//	p_n1_c_(p_n1_ / 2 + 1),
//	h_(_h),
//	off0_(1), off1_(1),
//	SCALE(1. / T(p_n0_ * p_n1_))
//{
//	using namespace Util;
//	using namespace Debug;
//	using namespace Detail;
//	//------------------
//	temp1_r.resize(p_n0_ * p_n1_);
//	temp2_r.resize(p_n0_ * p_n1_);
//	temp_c.resize(p_n0_ * p_n1_c_);
//	//------------------
//	sobel_x_c.resize(p_n0_ * p_n1_c_);
//	sobel_y_c.resize(p_n0_ * p_n1_c_);
//	//------------------
//	plan_r2c = FFTWPlanPtr<T>(FFTWPlaner<T>::r2c(p_n0_, p_n1_, temp1_r.data(), temp_c.data()));
//	plan_c2r = FFTWPlanPtr<T>(FFTWPlaner<T>::c2r(p_n0_, p_n1_, temp_c.data(), temp1_r.data()));
//	//---------------
//	std::fill(temp1_r.begin(), temp1_r.end(), 0);
//	temp1_r[0 + 0 * p_n1_] = -1.;
//	temp1_r[1 + 0 * p_n1_] = 0.;
//	temp1_r[2 + 0 * p_n1_] = 1.;
//	temp1_r[0 + 1 * p_n1_] = -2.;
//	temp1_r[1 + 1 * p_n1_] = 0.;
//	temp1_r[2 + 1 * p_n1_] = 2.;
//	temp1_r[0 + 2 * p_n1_] = -1.;
//	temp1_r[1 + 2 * p_n1_] = 0.;
//	temp1_r[2 + 2 * p_n1_] = 1.;
//	FFTWExecutor<T>::r2c(plan_r2c.get(), temp1_r.data(), sobel_x_c.data());
//	//---------------
//	temp1_r[0 + 0 * p_n1_] = -1.;
//	temp1_r[1 + 0 * p_n1_] = -2.;
//	temp1_r[2 + 0 * p_n1_] = -1.;
//	temp1_r[0 + 1 * p_n1_] = 0.;
//	temp1_r[1 + 1 * p_n1_] = 0.;
//	temp1_r[2 + 1 * p_n1_] = 0.;
//	temp1_r[0 + 2 * p_n1_] = 1.;
//	temp1_r[1 + 2 * p_n1_] = 2.;
//	temp1_r[2 + 2 * p_n1_] = 1.;
//	FFTWExecutor<T>::r2c(plan_r2c.get(), temp1_r.data(), sobel_y_c.data());
//}//TurboPacker::Spectral::Impl::Config::Config
//
//template<class T_, bool SIMD_, bool DEBUG_>
//std::ostream& TurboPacker::Spectral::Impl::Config<T_, SIMD_, DEBUG_>::operator<<(std::ostream& s) {
//	s << "----------- Config -----------" << std::endl;
//	s << "n0: " << n0_ << ", n1: " << n1_ << std::endl;
//	s << "p_n0: " << p_n0_ << ", p_n1_c: " << p_n1_c_ << std::endl;
//	s << "off0: " << off0_ << ", off1: " << off1_ << std::endl;
//	s << "Size: " << ((sizeof(std::complex<T>) * (sobel_x_c.size() + sobel_y_c.size() + 
//		temp_c.size()) + sizeof(T) * (temp1_r.size() + temp2_r.size())) / 1000000) << "mb" << std::endl;
//	s << "------------------------------" << std::endl;
//}//TurboPacker::Spectral::Impl::Config::operator<<
//
//template<class T_, bool SIMD_, bool DEBUG_>
//void TurboPacker::Spectral::Impl::fft_signal(
//	Config<T_, SIMD_, DEBUG_>& _c,
//	const FFTWVector<T_>& _map,
//	const FFTWVector<T_>& _map_log,
//	FFTWVector<std::complex<T_>>& _out,
//	FFTWVector<std::complex<T_>>& _out_log
//) {
//	using namespace Util;
//	using namespace Debug;
//	using namespace Detail;
//
//	if constexpr (DEBUG_) image_real<T_, false>(_c.p_n0_, _c.p_n1_, const_cast<T_*>(_map.data()), "1_signal_linear.png");
//	if constexpr (DEBUG_) image_real<T_, false>(_c.p_n0_, _c.p_n1_, const_cast<T_*>(_map_log.data()), "1_signal_log.png");
//
//	FFTWExecutor<T_>::r2c(_c.plan_r2c.get(), const_cast<T_*>(_map.data()), _out.data());
//	FFTWExecutor<T_>::r2c(_c.plan_r2c.get(), const_cast<T_*>(_map_log.data()), _out_log.data());
//
//}//TurboPacker::Spectral::Impl::correlate_signal
//
//template<class T_, bool SIMD_, bool DEBUG_>
//void TurboPacker::Spectral::Impl::gradients_signal(
//	Config<T_, SIMD_, DEBUG_>& _c,
//	const FFTWVector<std::complex<T_>>& _in,
//	FFTWVector<std::complex<T_>>& _out
//) {
//	using namespace Util;
//	using namespace Debug;
//	using namespace Detail;
//
//	if constexpr (!SIMD_) {
//
//		//------------------ Sobel X --------------------
//		for (int32 n0 = 0; n0 < _c.p_n0_; ++n0) {
//			for (int32 n1 = 0; n1 < _c.p_n1_c_; ++n1) {
//				const int32 i = n1 + n0 * _c.p_n1_c_;
//				const std::complex<T_>& c1 = _in[i];
//				const std::complex<T_>& c2 = _c.sobel_x_c[i];
//				const auto r = c1 * c2;
//				_c.temp_c[i] = r;
//			}
//		}
//		FFTWExecutor<T_>::c2r(_c.plan_c2r.get(), _c.temp_c.data(), _c.temp1_r.data());
//		for (int32 n0 = 0; n0 < _c.n0_; ++n0) {
//			for (int32 n1 = 0; n1 < _c.n1_; ++n1) {
//				const int32 i1 = (n1 + _c.off1_) + (n0 + _c.off0_) * _c.p_n1_;
//				const int32 i2 = (n1 + _c.off1_ + 1) + (n0 + _c.off0_ + 1) * _c.p_n1_; // +1 for kernel offset
//				const int32 v = int32(std::abs((_c.temp1_r[i2] * _c.SCALE)));
//				_c.temp2_r[i1] = v == 0 ? 0. : 1.;
//			}
//		}
//
//		//------------------ Sobel Y --------------------
//		for (int32 n0 = 0; n0 < _c.p_n0_; ++n0) {
//			for (int32 n1 = 0; n1 < _c.p_n1_c_; ++n1) {
//				const int32 i = n1 + n0 * _c.p_n1_c_;
//				const std::complex<T_>& c1 = _in[i];
//				const std::complex<T_>& c2 = _c.sobel_y_c[i];
//				const auto r = c1 * c2;
//				_c.temp_c[i] = r;
//			}
//		}
//		FFTWExecutor<T_>::c2r(_c.plan_c2r.get(), _c.temp_c.data(), _c.temp1_r.data());
//		for (int32 n0 = 0; n0 < _c.n0_; ++n0) {
//			for (int32 n1 = 0; n1 < _c.n1_; ++n1) {
//				const int32 i1 = (n1 + _c.off1_) + (n0 + _c.off0_) * _c.p_n1_;
//				const int32 i2 = (n1 + _c.off1_ + 1) + (n0 + _c.off0_ + 1) * _c.p_n1_;
//				const int32 v1 = int32(std::abs((_c.temp1_r[i2] * _c.SCALE)));
//				const int32 v2 = int32(_c.temp2_r[i1]);
//				_c.temp2_r[i1] = v2 == 0 && v1 == 0 ? 0. : 1.;
//			}
//		}
//
//		if constexpr (DEBUG_) image_real<T_, false>(_c.p_n0_, _c.p_n1_, _c.temp2_r.data(), "2_signal_gradients.png");
//
//		FFTWExecutor<T_>::r2c(_c.plan_r2c.get(), _c.temp2_r.data(), _out.data());
//	}
//
//}//TurboPacker::Spectral::Impl::gradients_signal
//
//template<class T_, class VOXLER_, bool SIMD_, bool DEBUG_>
//void TurboPacker::Spectral::Impl::fft_kernel(
//	Config<T_, SIMD_, DEBUG_>& _c,
//	const VOXLER_& _v,
//	FFTWVector<std::complex<T_>>& _out
//) {
//	using namespace Util;
//	using namespace Debug;
//	using namespace Detail;
//
//	if constexpr (!SIMD_) {
//		_v.fill(_c.temp1_r);
//		if constexpr (DEBUG_) image_real<T_, false>(_c.p_n0_, _c.p_n1_, _c.temp1_r.data(), "3_kernel_voxel.png");
//		FFTWExecutor<T_>::r2c(_c.plan_r2c.get(), _c.temp1_r.data(), _out.data());
//	}
//}//TurboPacker::Spectral::Impl::fft_kernel
//
//template<class T_, bool SIMD_, bool DEBUG_>
//void TurboPacker::Spectral::Impl::correlate_kernel_signal(
//	Config<T_, SIMD_, DEBUG_>& _c,
//	const FFTWVector<std::complex<T_>>& _kernel_in,
//	const FFTWVector<std::complex<T_>>& _signal_in,
//	const FFTWVector<std::complex<T_>>& _signal_log_in,
//	FFTWVector<T_>& _out,
//	FFTWVector<T_>& out_l
//) {
//	using namespace Util;
//	using namespace Debug;
//	using namespace Detail;
//
//	if constexpr (!SIMD_) {
//		//const int32 wr = int32(std::floor(_ext[0])) + 1;
//		//const int32 hr = int32(std::floor(_ext[1])) + 1;
//		//const int32 w = wr % 2 == 0 ? wr + 1 : wr;
//		//const int32 h = hr % 2 == 0 ? hr + 1 : hr;
//		//const int32 zh = int32(std::floor(_ext[2])) + 1;
//		//const int32 w2 = w / 2;
//		//const int32 h2 = h / 2;
//
//		////------------ Convolution with Signal -----------
//		//for (int32 n0 = 0; n0 < _c.p_n0_; ++n0) {
//		//	for (int32 n1 = 0; n1 < _c.p_n1_c_; ++n1) {
//		//		const int32 i = n1 + n0 * _c.p_n1_c_;
//		//		const std::complex<T_>& c1 = _signal_in[i];
//		//		const std::complex<T_>& c2 = _kernel_in[i];
//		//		const auto r = c1 * c2;
//		//		_c.temp_c[i] = r;
//		//	}
//		//}
//
//		//FFTWExecutor<T_>::c2r(_c.plan_c2r.get(), _c.temp_c.data(), _out.data());
//
//		//if constexpr (DEBUG_) image_real<T_, false>(_c.p_n0_, _c.p_n1_, _c.temp1_r.data(), "4_convolution.png");
//
//		////------------ Convolution with Signal -----------
//		//for (int32 n0 = 0; n0 < _c.n0_; ++n0) {
//		//	for (int32 n1 = 0; n1 < _c.n1_; ++n1) {
//		//		const int32 i1 = n1 + n0 * _c.n1_;
//		//		const int32 i2 = (n1 + _c.off1_ + h2) + (n0 + _c.off0_ + w2) * _c.p_n1_;
//		//		const int32 i3 = (n1 + _c.off1_) + (n0 + _c.off0_) * _c.p_n1_;
//		//		_out[i3] = temp1_r[i2] * SCALE;
//		//	}
//		//}
//		////------------ Log Convolution with Signal -----------
//		//for (int32 n0 = 0; n0 < _c.p_n0_; ++n0) {
//		//	for (int32 n1 = 0; n1 < _c.p_n1_c_; ++n1) {
//		//		const int32 i = n1 + n0 * _c.p_n1_c_;
//		//		const std::complex<T_>& c1 = _signal_log_in[i];
//		//		const std::complex<T_>& c2 = _kernel_in[i];
//		//		const auto r = c1 * c2;
//		//		_c.temp_c[i] = r;
//		//	}
//		//}
//		//FFTWExecutor<T>::c2r(_c.plan_c2r.get(), _c.temp_c.data(), out_l.data());
//
//		////--------- Rectified solution --------
//		//for (int32 n0 = 0; n0 < _c.n0_; ++n0) {
//		//	for (int32 n1 = 0; n1 < _c.n1_; ++n1) {
//		//		const int32 i1 = n1 + n0 * _c.n1_;
//		//		const int32 i2 = (n1 + _c.off1_ + h2) + (n0 + _c.off0_ + w2) * _c.p_n1_;
//		//		const int32 i3 = (n1 + _c.off1_) + (n0 + _c.off0_) * _c.p_n1_;
//		//		out_l[i3] = temp1_r[i2] * _c.SCALE;
//		//	}
//		//}
//		//if constexpr (DEBUG_) image_real<T_, false>(_c.p_n0_, _c.p_n1_, _c.temp2_r.data(), "5_solution.png");
//	}
//}//TurboPacker::Spectral::Impl::correlate_kernel_signal
//
//template<class T_, bool SIMD_, bool DEBUG_>
//void TurboPacker::Spectral::Impl::extract_solution(
//	Config<T_, SIMD_, DEBUG_>& _c,
//	const FFTWVector<T_>& _corr_in,
//	const FFTWVector<T_>& _corr_log_in,
//	FFTWVector<T_>& _sol_out
//) {
//	using namespace Util;
//	using namespace Debug;
//	using namespace Detail;
//
//	if constexpr (!SIMD_) {
//
//		////--------- Rectified solution --------
//		for (int32 n0 = 0; n0 < n0_; ++n0) {
//			for (int32 n1 = 0; n1 < n1_; ++n1) {
//				const int32 i1 = n1 + n0 * n1_;
//				const int32 i2 = (n1 + off1_ + h2) + (n0 + off0_ + w2) * p_n1_;
//				const int32 i3 = (n1 + off1_) + (n0 + off0_) * p_n1_;
//
//				const T expec = T(w * h) * std::log(_c.map[i1] + 1.);
//				const T r = std::abs((_c.temp1_r[i2] * SCALE));
//				_c.temp2_r[i3] = (FMath::IsNearlyEqual(temp2_r[i3], 1.) && FMath::IsNearlyEqual(expec, r)) ? 1. : 0.;
//			}
//		}
//		if constexpr (DEBUG) image_real<T, false>(p_n0_, p_n1_, temp2_r.data(), "5_solution.png");
//
//		
//	}
//}//TurboPacker::Spectral::Impl::extract_solution
//
//template<class T_, bool SIMD_, bool DEBUG_>
//void TurboPacker::Spectral::Impl::gradients_signal_stencil(
//	Config<T_, SIMD_, DEBUG_>& _config,
//	const FFTWVector<T_>& _sol,
//	FFTWVector<std::complex<T_>>& _border_out
//) {
//	using namespace Util;
//	using namespace Debug;
//	using namespace Detail;
//
//	if constexpr (!SIMD_) {
//		for (int32 n0 = 0; n0 < n0_; ++n0) {
//			for (int32 n1 = 0; n1 < n1_; ++n1) {
//
//				const int32 i0 = n1 + n0 * _c.n1_;
//
//				const int32 i1 = (n1 - 1) + n0 * _c.n1_;
//				const int32 i2 = n1 + (n0 + 1) * _c.n1_;
//				const int32 i3 = (n1 + 1) + n0 * _c.n1_;
//				const int32 i4 = n1 + (n0 - 1) * _c.n1_;
//
//				_config.temp1_r[i0] = int32(_sol[i1]) + int32(_sol[i3]) == 1 ||
//					int32(_sol[i2]) + int32(_sol[i4]) == 1 ? 1. : 0.;
//
//			}
//		}
//
//		if constexpr (DEBUG_) image_real<T_, false>(_c.p_n0_, _c.p_n1_, _c.temp1_r.data(), "6_border.png");
//		FFTWExecutor<T_>::r2c(_c.plan_r2c.get(), _c.temp1_r.data(), _border_out.data());
//
//	}
//
//}//TurboPacker::Spectral::Impl::gradients_signal_stencil
//
////------------------------------------------------
////------------------ MedianQuadTree --------------
////------------------------------------------------
//// needs debugging!
//
//template<class T, bool SIMD, bool DEBUG>
//void TurboPacker::Spectral::MedianQuadTree<T, SIMD, DEBUG>::Bucket::recompute() {
//	// n1: x, n0: y
//	max = -std::numeric_limits<T>::infinity();
//	min = std::numeric_limits<T>::infinity();
//
//	std::vector<std::pair<T, FIntVector2>> m;
//	m.reserve((bmax.X - bmin.X + 1) * (bmax.Y - bmin.Y + 1));
//
//	for (int32 n0 = bmin.Y; n0 < bmax.Y; ++n0) {
//		for (int32 n1 = bmin.X; n1 < bmax.X; ++n1) {
//			const int32 i = n1 + n0 * p->n1_;
//			if (i > p->map.size()) continue;
//			m.push_back({ p->map[i], { n0, n1} });
//			max = std::max(max, p->map[i]);
//			min = std::min(max, p->map[i]);
//		}
//	}
//
//	std::sort(m.begin(), m.end(), [](const auto& _e1, const auto& _e2) {
//		return std::get<0>(_e1) < std::get<0>(_e2);
//	});
//
//	if (m.size() % 2 == 0) {
//		const int32 idx = (m.size() + 1) / 2;
//		median = std::get<0>(m[idx]);
//
//		l.clear();
//		for (int32 j = 0; j <= idx; ++j)
//			l.push_back(std::get<1>(m[j]));
//
//		h.clear();
//		for (int32 j = idx + 1; j < m.size(); ++j)
//			h.push_back(std::get<1>(m[j]));
//
//	} else {
//		const int32 idx1 = m.size() / 2;
//		const int32 idx2 = idx1 + 1;
//		median = (std::get<0>(m[idx1]) + std::get<0>(m[idx2])) * 0.5;
//
//		l.clear();
//		for (int32 j = 0; j <= idx1; ++j)
//			l.push_back(std::get<1>(m[j]));
//
//		h.clear();
//		for (int32 j = idx2; j < m.size(); ++j)
//			h.push_back(std::get<1>(m[j]));
//	}
//
//	if constexpr (DEBUG) {
//		std::cout << level << " - " << median << std::endl;
//		std::cout << l.size() << " - " << h.size() << std::endl;
//	}
//
//}//Bucket::recompute
//
//template<class T, bool SIMD, bool DEBUG>
//std::pair<int32, int32> TurboPacker::Spectral::MedianQuadTree<T, SIMD, DEBUG>::Bucket::overlap(
//	const FIntVector2& _min, 
//	const FIntVector2& _max, 
//	T _h
//) const {
//
//	const auto contains = [&](int32 _x, int32 _y) {
//		return _min.X <= _x && _x < _max.X && _min.Y <= _y && _y < _max.Y;
//	};
//
//	//partial overlap
//	if (_min.X > bmin.X || _max.X < bmax.X || _min.Y > bmin.Y || _max.Y < bmax.Y) {
//
//		int32 low = 0;
//		int32 high = 0;
//
//		for (int32 n0 = std::max(bmin.Y, _min.Y); n0 < std::min(bmax.Y, _max.Y); ++n0) {
//			for (int32 n1 = std::max(bmin.X, _min.X); n1 < std::min(bmax.X, _max.X); ++n1) {
//				const int32 i = n0 + n1 * p->n1_;
//				if (i > p->map.size()) continue;
//				if (p->map[i] >= _h) high++;
//				else low++;
//			}
//		}
//
//		return { low, high };
//
//	}
//
//	//--------------------
//
//	if (isFlat) {
//
//		if (_h >= median) return { 0, l.size() + h.size() };
//		else return { l.size() + h.size(), 0 };
//
//	}
//	
//	//--------------------
//	
//	if (_h > median) {
//
//		int32 low = l.size();
//		int32 high = 0;
//
//		for (const FIntVector2 pos : h) {
//			const int32 i = pos.Y + pos.X * p->n1_;
//			if (!contains(pos.X, pos.Y)) continue;
//
//			if (p->map[i] >= _h)
//				high++;
//			if (p->map[i] < _h)
//				low++;
//		}
//
//		return { low, high };
//
//	} else {
//
//		int32 low = 0;
//		int32 high = h.size();
//
//		for (const FIntVector2 pos : l) {
//			const int32 i = pos.Y + pos.X * p->n1_;
//			if (!contains(pos.X, pos.Y)) continue;
//
//			if (p->map[i] >= _h)
//				high++;
//			if (p->map[i] < _h)
//				low++;
//		}
//
//		return { low, high };
//
//	}
//
//}//Bucket::overlap
//
////--------------
//
//template<class T, bool SIMD, bool DEBUG>
//void TurboPacker::Spectral::MedianQuadTree<T, SIMD, DEBUG>::Node::recompute() {
//
//	for (int32 i = 0; i < 4; ++i){
//		switch (c[0].index()) {
//			case 0:
//				std::get<0>(c[i])->recompute();
//			break;
//			case 1:
//				std::get<1>(c[i])->recompute();
//			break;
//		}
//	}
//
//	//--------------
//
//	std::vector<std::pair<int32, T>> m;
//
//	switch (c[0].index()) {
//		case 0:
//		{			
//			for (int32 i = 0; i < 4; ++i)
//				m.push_back({ i, std::get<0>(c[i])->median });
//		}
//		break;
//		case 1:
//		{
//			for (int32 i = 0; i < 4; ++i)
//				m.push_back({ i, std::get<1>(c[i])->median });
//		}
//		break;
//	}
//
//	std::sort(m.begin(), m.end(), [](const auto& _e1, const auto& _e2) {
//		return std::get<1>(_e1) < std::get<1>(_e2);
//	});
//
//	median = (std::get<1>(m[1]) + std::get<1>(m[2])) * 0.5;
//	min = std::get<1>(m[0]);
//	max = std::get<1>(m[3]);
//
//	l[0] = std::get<0>(m[0]);
//	l[1] = std::get<0>(m[1]);
//
//	h[0] = std::get<0>(m[2]);
//	h[1] = std::get<0>(m[3]);
//
//	if constexpr (DEBUG) {
//		std::cout << level << " - " << median << std::endl;
//	}
//
//}//Node::recompute
//
//template<class T, bool SIMD, bool DEBUG>
//std::pair<int32, int32> TurboPacker::Spectral::MedianQuadTree<T, SIMD, DEBUG>::Node::overlap(
//	const FIntVector2& _min,
//	const FIntVector2& _max,
//	T _h
//) const {
//
//	//partial overlap
//	if (_min.X > bmin.X || _max.X < bmax.X || _min.Y > bmin.Y || _max.Y < bmax.Y) {
//		int32 low = 0;
//		int32 high = 0;
//		switch (c[0].index()) {
//			case 0:
//			{
//				for (int32 i = 0; i < 4; ++i) {
//					const auto r = std::get<0>(c[i])->overlap(_min, _max, _h);
//					low += std::get<0>(r);
//					high += std::get<1>(r);
//				}
//			}
//			break;
//			case 1:
//			{
//				for (int32 i = 0; i < 4; ++i) {
//					const auto r = std::get<1>(c[i])->overlap(_min, _max, _h);
//					low += std::get<0>(r);
//					high += std::get<1>(r);
//					//std::cout << std::get<0>(r) << ", " << std::get<1>(r) << std::endl;
//				}
//			}
//			break;
//		}
//		return { low, high };
//	}
//
//	//--------------------
//
//	if (isFlat) {
//
//		if (_h >= median) return { 0, 4 * (p->maxLevel - level) * p->extend * p->extend };
//		else return { 4 * (p->maxLevel - level) * p->extend * p->extend, 0 };
//
//	}
//
//	//--------------------
//	
//	if (_h > median) {
//
//		int32 low = l.size();
//		int32 high = 0;
//
//		switch (c[0].index()) {
//			case 0:
//			{
//				const auto r1 = std::get<0>(c[h[0]])->overlap(_min, _max, _h);
//				const auto r2 = std::get<0>(c[h[1]])->overlap(_min, _max, _h);
//				low += std::get<0>(r1) + std::get<0>(r2);
//				high += std::get<1>(r1) + std::get<1>(r2);
//			}
//			break;
//			case 1:
//			{
//				const auto r1 = std::get<1>(c[h[0]])->overlap(_min, _max, _h);
//				const auto r2 = std::get<1>(c[h[1]])->overlap(_min, _max, _h);
//				low += std::get<0>(r1) + std::get<0>(r2);
//				high += std::get<1>(r1) + std::get<1>(r2);
//			}
//			break;
//		}
//		
//		return { low, high };
//
//	} else {
//
//		int32 low = 0;
//		int32 high = h.size();
//
//		switch (c[0].index()) {
//			case 0:
//			{
//				const auto r1 = std::get<0>(c[l[0]])->overlap(_min, _max, _h);
//				const auto r2 = std::get<0>(c[l[1]])->overlap(_min, _max, _h);
//				low += std::get<0>(r1) + std::get<0>(r2);
//				high += std::get<1>(r1) + std::get<1>(r2);
//			}
//			break;
//			case 1:
//			{
//				const auto r1 = std::get<1>(c[l[0]])->overlap(_min, _max, _h);
//				const auto r2 = std::get<1>(c[l[1]])->overlap(_min, _max, _h);
//				low += std::get<0>(r1) + std::get<0>(r2);
//				high += std::get<1>(r1) + std::get<1>(r2);
//			}
//			break;
//		}
//
//		return { low, high };
//
//	}
//
//}//Node::overlap
//
////--------------
//
//template<class T, bool SIMD, bool DEBUG>
//TurboPacker::Spectral::MedianQuadTree<T, SIMD, DEBUG>::MedianQuadTree(
//	Impl::Config<T, SIMD, DEBUG>& _c,
//	const FFTWVector<T>& _map,
//	const int32 _min_ext
//) : map(_map), n0_(_c.n0_), n1_(_c.n1_), extend(_min_ext), world(_c.world) {
//
//	const int32 s = std::max(_c.n0_, _c.n1_);
//	maxLevel = 1;
//	int32 t = _min_ext;
//	while (t < s) {
//		maxLevel++;
//		t *= 2;
//	}
//	if constexpr (DEBUG) {
//		std::cout << "s: " << s << std::endl;
//		std::cout << "t: " << t << std::endl;
//		std::cout << "maxLevel: " << maxLevel << std::endl;
//		std::cout << "----------------" << std::endl;
//	}
//	//--------------
//	root = std::make_unique<Node>(this, FIntVector2(0), FIntVector2(t), 1);
//
//	std::queue<Node*> q;
//	q.push(root.get());
//
//	int32 nc = 0;
//	int32 bc = 0;
//	while (!q.empty()) {
//
//		Node* n = q.front();
//		q.pop();
//
//		nc++;
//
//		/*
//		children:
//			0 | 1
//			2 | 3
//		*/
//		const FIntVector2 min_0 = FIntVector2(n->bmin.X, n->bmin.Y + (n->bmax.Y - n->bmin.Y) / 2);
//		const FIntVector2 max_0 = FIntVector2(n->bmin.X + (n->bmax.X - n->bmin.X) / 2, n->bmax.Y);
//
//		const FIntVector2 min_1 = FIntVector2(n->bmin.X + (n->bmax.X - n->bmin.X) / 2, n->bmin.Y + (n->bmax.Y - n->bmin.Y) / 2);
//		const FIntVector2 max_1 = n->bmax;
//
//		const FIntVector2 min_2 = n->bmin;
//		const FIntVector2 max_2 = FIntVector2(n->bmin.X + (n->bmax.X - n->bmin.X) / 2, n->bmin.Y + (n->bmax.Y - n->bmin.Y) / 2);
//
//		const FIntVector2 min_3 = FIntVector2(n->bmin.X + (n->bmax.X - n->bmin.X) / 2, n->bmin.Y);
//		const FIntVector2 max_3 = FIntVector2(n->bmax.X, n->bmin.Y + (n->bmax.Y - n->bmin.Y) / 2);
//
//		//std::cout << "---- " << n->level + 1 << " ----" << std::endl;
//		//std::cout << n->bmin << n->bmax << std::endl;
//		//std::cout << "0: " << min_0 << max_0 << std::endl;
//		//std::cout << "1: " << min_1 << max_1 << std::endl;
//		//std::cout << "2: " << min_2 << max_2 << std::endl;
//		//std::cout << "3: " << min_3 << max_3 << std::endl;
//		//std::cout << "----------------" << std::endl;
//
//		if constexpr (DEBUG) {
//			DrawDebugBox(world,
//				FVector(n->bmin.X + (n->bmax.X - n->bmin.X) / 2, n->bmin.Y + (n->bmax.Y - n->bmin.Y) / 2, (n->bmax.X - n->bmin.X) / 2),
//				FVector((n->bmax.X - n->bmin.X) / 2),
//				FColor::Blue, true);
//		}
//
//		if (n->level == maxLevel-1) {
//
//			bc += 4;
//
//			// --- 0 ---
//			n->c[0] = std::make_unique<Bucket>(this, min_0, max_0, n->level + 1);
//			
//			// --- 1 ---
//			n->c[1] = std::make_unique<Bucket>(this, min_1, max_1, n->level + 1);
//
//			// --- 2 ---
//			n->c[2] = std::make_unique<Bucket>(this, min_2, max_2, n->level + 1);
//
//			// --- 3 ---
//			n->c[3] = std::make_unique<Bucket>(this, min_3, max_3, n->level + 1);
//
//			if constexpr (DEBUG) {
//
//				DrawDebugBox(world,
//					FVector(min_0.X + (max_0.X - min_0.X) / 2, min_0.Y + (max_0.Y - min_0.Y) / 2, (max_0.X - min_0.X) / 2),
//					FVector((max_0.X - min_0.X) / 2),
//					FColor::Green, true);
//
//				DrawDebugBox(world,
//					FVector(min_1.X + (max_1.X - min_1.X) / 2, min_1.Y + (max_1.Y - min_1.Y) / 2, (max_1.X - min_1.X) / 2),
//					FVector((max_1.X - min_1.X) / 2),
//					FColor::Green, true);
//
//				DrawDebugBox(world,
//					FVector(min_2.X + (max_2.X - min_2.X) / 2, min_2.Y + (max_2.Y - min_2.Y) / 2, (max_2.X - min_2.X) / 2),
//					FVector((max_2.X - min_2.X) / 2),
//					FColor::Green, true);
//
//				DrawDebugBox(world,
//					FVector(min_3.X + (max_3.X - min_3.X) / 2, min_3.Y + (max_3.Y - min_3.Y) / 2, (max_3.X - min_3.X) / 2),
//					FVector((max_3.X - min_3.X) / 2),
//					FColor::Green, true);
//			}
//
//		} else {
//
//			n->c[0] = std::make_unique<Node>(this, min_0, max_0, n->level + 1);
//			q.push(std::get<0>(n->c[0]).get());
//
//			n->c[1] = std::make_unique<Node>(this, min_1, max_1, n->level + 1);
//			q.push(std::get<0>(n->c[1]).get());
//
//			n->c[2] = std::make_unique<Node>(this, min_2, max_2, n->level + 1);
//			q.push(std::get<0>(n->c[2]).get());
//
//			n->c[3] = std::make_unique<Node>(this, min_3, max_3, n->level + 1);
//			q.push(std::get<0>(n->c[3]).get());
//
//		}
//	}
//
//	if constexpr (DEBUG) {
//		std::cout << "Nodes: " << nc << std::endl;
//		std::cout << "Buckets: " << bc << std::endl;
//	}
//
//}//TurboPacker::Spectral::MedianQuadTree::MedianQuadTree
//
//template<class T, bool SIMD, bool DEBUG>
//void TurboPacker::Spectral::MedianQuadTree<T, SIMD, DEBUG>::recompute() {
//	root->recompute();
//}//TurboPacker::Spectral::MedianQuadTree::push
//
//template<class T, bool SIMD, bool DEBUG>
//std::pair<int32, int32> TurboPacker::Spectral::MedianQuadTree<T, SIMD, DEBUG>::check_height(
//	const FIntVector& _pos, 
//	const FIntVector2& _ext
//) const {
//	if constexpr (DEBUG) {
//		DrawDebugBox(world,
//			FVector(_pos.X, _pos.Y, _pos.Z + 5.),
//			FVector(_ext.X, _ext.Y, 5.),
//			FColor::Red, true);
//	}
//	return root->overlap(
//		FIntVector2(_pos.X, _pos.Y) - _ext,
//		FIntVector2(_pos.X, _pos.Y) + _ext,
//		_pos.Z);
//}//TurboPacker::Spectral::MedianQuadTree::check_height
//
////------------------------------------------------
////------------------ BoxStack --------------------
////------------------------------------------------
//
//template<class T>
//void TurboPacker::Spectral::BoxStack<T>::push(const UE::Math::TBox<T>& _box) {
//
//	const int32 min = int32(std::round(_box.Min.Z));
//	const int32 max = int32(std::round(_box.Max.Z));
//
//	boxes.emplace(++idd, _box);
//	list.emplace_back(min, max, idd);
//
//	boxs.push_back(_box);
//
//}//TurboPacker::Spectral::BoxStack::push
//
//template<class T>
//bool TurboPacker::Spectral::BoxStack<T>::overlap(const UE::Math::TBox<T>& _box) const {
//
//	const int32 min = int32(std::round(_box.Min.Z));
//	const int32 max = int32(std::round(_box.Max.Z));
//
//	for (const auto& [h_min, h_max, id] : list) {
//		if (max < h_min || h_max < min) {
//			const auto& b = boxes.find(id).value();
//			if (_box.Intersect(b)) return true;
//		}
//	}
//
//	return false;
//}//TurboPacker::Spectral::BoxStack::overlap
//
//template<class T>
//bool TurboPacker::Spectral::BoxStack<T>::overlap_naive(const UE::Math::TBox<T>& _box) const {
//
//	for (const auto b : boxs)
//		if (_box.Intersect(b)) return true;
//
//	return false;
//}//TurboPacker::Spectral::BoxStack::overlap
//
////------------------------------------------------
////------------------ HeightMap -------------------
////------------------------------------------------
//
//template<class T, bool SIMD, bool DEBUG>
//TurboPacker::Spectral::HeightMap<T, SIMD, DEBUG>::HeightMap(const int32 _n0, const int32 _n1, const int32 _h) {
//	using namespace Util;
//	using namespace Detail;
//	//-------------------------
//	config_ = std::make_unique<Impl::Config<T, SIMD, DEBUG>>(_n0, _n1, _h);
//	//-------------------------
//	map.resize(config_->p_n0_ * config_->p_n1_);
//	map_log.resize(config_->p_n0_ * config_->p_n1_);
//	std::fill(map.begin(), map.end(), 0);
//	std::fill(map_log.begin(), map_log.end(), 1.);
//	//-------------------------
//	//border n0
//	for (int32 n0 = 0; n0 <= config_->n0_ + 1; ++n0) {
//		const int32 i1 = n0 * config_->p_n1_;
//		const int32 i2 = config_->n1_ + 1 + n0 * config_->p_n1_;
//		map[i1] = T(config_->h_);
//		map[i2] = T(config_->h_);
//		map_log[i1] = std::log(T(config_->h_) + 1);
//		map_log[i2] = std::log(T(config_->h_) + 1);
//	}
//	//border n1
//	for (int32 n1 = 0; n1 <= config_->n1_ + 1; ++n1) {
//		const int32 i1 = n1;
//		const int32 i2 = n1 + (config_->n0_ + 1) * config_->p_n1_;
//		map[i1] = T(config_->h_);
//		map[i2] = T(config_->h_);
//		map_log[i1] = std::log(T(config_->h_) + 1);
//		map_log[i2] = std::log(T(config_->h_) + 1);
//	}
//	//-------------------------
//	signal_c.resize(config_->p_n0_ * config_->p_n1_c_);
//	signal_log_c.resize(config_->p_n0_ * config_->p_n1_c_);
//	kernel_c.resize(config_->p_n0_ * config_->p_n1_c_);
//	border_c.resize(config_->p_n0_ * config_->p_n1_c_);
//	//-------------------------
//	costf = [](int32 _n0, int32 _n1, int32 _h, int32 _c) {
//		return std::pow<T>(T(_h), 3) + T(_n0) + T(_n1) - _c;
//	};
//}//TurboPacker::Spectral::HeightMap::HeightMap
//
//template<class T, bool SIMD, bool DEBUG>
//std::vector<std::tuple<FIntVector, T, int32>> TurboPacker::Spectral::HeightMap<T, SIMD, DEBUG>::overlap(const FVector& _ext, const int32 _perm) {
//	using namespace Util;
//	using namespace Debug;
//	using namespace Detail;
//	using namespace Impl;
//
//	fft_signal<T, SIMD, DEBUG>(*config_, map, map_log, signal_c, signal_log_c);
//
//
//	FMath::IsNearlyEqual();
//
//
//
//
//
//
//
//
//
//
//
//
//	//------------------------------------------------
//	//-------------- Signal/Convolution --------------
//	//------------------------------------------------
//
//	//for (int32 n0 = 0; n0 < n0_; ++n0) {
//	//	for (int32 n1 = 0; n1 < n1_; ++n1) {	
//	//		const int32 i1 = n1 + n0 * n1_;
//	//		const int32 i2 = (n1 + off1_) + (n0 + off0_) * p_n1_;
//	//		temp1_r[i2] = T(map[i1]);
//	//	}
//	//}
//	//for (int32 n0 = -1; n0 <= n0_; ++n0) {
//	//	const int32 i1 = (-1 + off1_) + (n0 + off0_) * p_n1_;
//	//	const int32 i2 = (n1_ + off1_) + (n0 + off0_) * p_n1_;
//	//	temp1_r[i1] = T(h_);
//	//	temp1_r[i2] = T(h_);
//	//}
//	//for (int32 n1 = -1; n1 <= n1_; ++n1) {
//	//	const int32 i1 = (n1 + off1_) + (-1 + off0_) * p_n1_;
//	//	const int32 i2 = (n1 + off1_) + (n0_ + off0_) * p_n1_;
//	//	temp1_r[i1] = T(h_);
//	//	temp1_r[i2] = T(h_);
//	//}
//	//if constexpr (DEBUG) image_real<T, false>(p_n0_, p_n1_, temp1_r.data(), "1_signal.png");
//	//FFTWExecutor<T>::r2c(plan_r2c.get(), temp1_r.data(), signal_c.data());
//
//	//for (int32 n0 = 0; n0 < n0_; ++n0) {
//	//	for (int32 n1 = 0; n1 < n1_; ++n1) {
//	//		const int32 i1 = n1 + n0 * n1_;
//	//		const int32 i2 = (n1 + off1_) + (n0 + off0_) * p_n1_;
//	//		temp1_r[i2] = std::log(T(map[i1]+1));
//	//	}
//	//}
//	//for (int32 n0 = -1; n0 <= n0_; ++n0) {
//	//	const int32 i1 = (-1 + off1_) + (n0 + off0_) * p_n1_;
//	//	const int32 i2 = (n1_ + off1_) + (n0 + off0_) * p_n1_;
//	//	temp1_r[i1] = std::log(T(h_)+1);
//	//	temp1_r[i2] = std::log(T(h_) + 1);
//	//}
//	//for (int32 n1 = -1; n1 <= n1_; ++n1) {
//	//	const int32 i1 = (n1 + off1_) + (-1 + off0_) * p_n1_;
//	//	const int32 i2 = (n1 + off1_) + (n0_ + off0_) * p_n1_;
//	//	temp1_r[i1] = std::log(T(h_) + 1);
//	//	temp1_r[i2] = std::log(T(h_) + 1);
//	//}
//	//FFTWExecutor<T>::r2c(plan_r2c.get(), temp1_r.data(), signal_log_c.data());
//
//	////------------------------------------------------
//	////-------------- Sobel/Convolution ---------------
//	////------------------------------------------------
//
//	////------------------ Sobel X --------------------
//	//for (int32 n0 = 0; n0 < p_n0_; ++n0) {
//	//	for (int32 n1 = 0; n1 < p_n1_c_; ++n1) {
//	//		const int32 i = n1 + n0 * p_n1_c_;
//	//		const std::complex<T>& c1 = signal_c[i];
//	//		const std::complex<T>& c2 = sobel_x_c[i];
//	//		const auto r = c1 * c2;
//	//		temp_c[i] = r;
//	//	}
//	//}
//	//FFTWExecutor<T>::c2r(plan_c2r.get(), temp_c.data(), temp1_r.data());
//	//for (int32 n0 = 0; n0 < n0_; ++n0) {
//	//	for (int32 n1 = 0; n1 < n1_; ++n1) {
//	//		const int32 i1 = (n1 + off1_) + (n0 + off0_) * p_n1_;
//	//		const int32 i2 = (n1 + off1_ + 1) + (n0 + off0_ + 1) * p_n1_;
//	//		const int32 v = int32(std::abs((temp1_r[i2] * SCALE)));
//	//		temp2_r[i1] = v == 0 ? 0. : 1.;
//	//	}
//	//}
//
//	////------------------ Sobel Y --------------------
//	//for (int32 n0 = 0; n0 < p_n0_; ++n0) {
//	//	for (int32 n1 = 0; n1 < p_n1_c_; ++n1) {
//	//		const int32 i = n1 + n0 * p_n1_c_;
//	//		const std::complex<T>& c1 = signal_c[i];
//	//		const std::complex<T>& c2 = sobel_y_c[i];
//	//		const auto r = c1 * c2;
//	//		temp_c[i] = r;
//	//	}
//	//}
//	//FFTWExecutor<T>::c2r(plan_c2r.get(), temp_c.data(), temp1_r.data());
//	//for (int32 n0 = 0; n0 < n0_; ++n0) {
//	//	for (int32 n1 = 0; n1 < n1_; ++n1) {
//	//		const int32 i1 = (n1 + off1_) + (n0 + off0_) * p_n1_;
//	//		const int32 i2 = (n1 + off1_ + 1) + (n0 + off0_ + 1) * p_n1_;
//	//		const int32 v1 = int32(std::abs((temp1_r[i2] * SCALE)));
//	//		const int32 v2 = int32(temp2_r[i1]);
//	//		temp2_r[i1] = v2 == 0 && v1 == 0 ? 0. : 1.;
//	//	}
//	//}
//	//if constexpr (DEBUG) image_real<T, false>(p_n0_, p_n1_, temp2_r.data(), "2_sobel_bit.png");
//	//FFTWExecutor<T>::r2c(plan_r2c.get(), temp2_r.data(), border_c.data());
//
//	////------------------------------------------------
//	////----------------- Permutations -----------------
//	////------------------------------------------------
//	//const int32 wr = int32(std::floor(_ext[0])) + 1;
//	//const int32 hr = int32(std::floor(_ext[1])) + 1;
//	//const int32 w = wr % 2 == 0 ? wr + 1 : wr;
//	//const int32 h = hr % 2 == 0 ? hr + 1 : hr;
//	//const int32 zh = int32(std::floor(_ext[2])) + 1;
//	//const int32 w2 = w / 2;
//	//const int32 h2 = h / 2;
//
//	////std::cout << w << ", " << h << ", " << zh << std::endl;
//
//	//const auto perm = [&](T _x, T _y, int32 _per)->void {
//
//	//	//------------------------------------------------
//	//	//------------------- Kernel ---------------------
//	//	//------------------------------------------------
//	//	std::fill(temp1_r.begin(), temp1_r.end(), 0);
//	//	for (int32 n0 = 0; n0 < h; ++n0) {
//	//		for (int32 n1 = 0; n1 < w; ++n1) {
//	//			const int32 i = n1 + n0 * p_n1_;
//	//			temp1_r[i] = 1.f;
//	//		}
//	//	}
//	//	if constexpr (DEBUG) image_real<T, false>(p_n0_, p_n1_, temp1_r.data(), "3_kernel.png");
//	//	FFTWExecutor<T>::r2c(plan_r2c.get(), temp1_r.data(), kernel_c.data());
//	//	//------------ Convolution with Signal -----------
//	//	for (int32 n0 = 0; n0 < p_n0_; ++n0) {
//	//		for (int32 n1 = 0; n1 < p_n1_c_; ++n1) {
//	//			const int32 i = n1 + n0 * p_n1_c_;
//	//			const std::complex<T>& c1 = signal_c[i];
//	//			const std::complex<T>& c2 = kernel_c[i];
//	//			const auto r = c1 * c2;
//	//			temp_c[i] = r;
//	//		}
//	//	}
//	//	FFTWExecutor<T>::c2r(plan_c2r.get(), temp_c.data(), temp1_r.data());
//	//	if constexpr (DEBUG) image_real<T, false>(p_n0_, p_n1_, temp1_r.data(), "4_convolution.png");
//	//	//------------ Convolution with Signal -----------
//	//	for (int32 n0 = 0; n0 < n0_; ++n0) {
//	//		for (int32 n1 = 0; n1 < n1_; ++n1) {
//	//			const int32 i1 = n1 + n0 * n1_;
//	//			const int32 i2 = (n1 + off1_ + h2) + (n0 + off0_ + w2) * p_n1_;
//	//			const int32 i3 = (n1 + off1_) + (n0 + off0_) * p_n1_;
//
//	//			const int32 expec = w * h * map[i1];
//	//			const int32 r = int32(std::abs(std::round(temp1_r[i2] * SCALE)));
//	//			//signal_r[i2] = r > expec ? 2. : r < expec ? 0. : 1.;
//	//			temp2_r[i3] = r == expec ? 1 : 0;
//	//			//result[_idx][i1] = r > expec ? -1 : r < expec ? 1 : 0;
//	//		}
//	//	}
//	//	//------------ Log Convolution with Signal -----------
//	//	for (int32 n0 = 0; n0 < p_n0_; ++n0) {
//	//		for (int32 n1 = 0; n1 < p_n1_c_; ++n1) {
//	//			const int32 i = n1 + n0 * p_n1_c_;
//	//			const std::complex<T>& c1 = signal_log_c[i];
//	//			const std::complex<T>& c2 = kernel_c[i];
//	//			const auto r = c1 * c2;
//	//			temp_c[i] = r;
//	//		}
//	//	}
//	//	FFTWExecutor<T>::c2r(plan_c2r.get(), temp_c.data(), temp1_r.data());
//	//	//--------- Rectified solution --------
//	//	for (int32 n0 = 0; n0 < n0_; ++n0) {
//	//		for (int32 n1 = 0; n1 < n1_; ++n1) {
//	//			const int32 i1 = n1 + n0 * n1_;
//	//			const int32 i2 = (n1 + off1_ + h2) + (n0 + off0_ + w2) * p_n1_;
//	//			const int32 i3 = (n1 + off1_) + (n0 + off0_) * p_n1_;
//
//	//			const T expec = T(w * h) * std::log(T(map[i1] + 1));
//	//			const T r = std::abs((temp1_r[i2] * SCALE));
//	//			temp2_r[i3] = (FMath::IsNearlyEqual(temp2_r[i3], 1.) && FMath::IsNearlyEqual(expec, r)) ? 1. : 0.;
//	//		}
//	//	}
//	//	if constexpr (DEBUG) image_real<T, false>(p_n0_, p_n1_, temp2_r.data(), "5_solution.png");
//	//	//--------- Convolution of Kernel with border --------
//	//	for (int32 n0 = 0; n0 < p_n0_; ++n0) {
//	//		for (int32 n1 = 0; n1 < p_n1_c_; ++n1) {
//	//			const int32 i = n1 + n0 * p_n1_c_;
//	//			const std::complex<T>& c1 = border_c[i];
//	//			const std::complex<T>& c2 = kernel_c[i];
//	//			const auto r = c1 * c2;
//	//			temp_c[i] = r;
//	//		}
//	//	}
//	//	FFTWExecutor<T>::c2r(plan_c2r.get(), temp_c.data(), temp1_r.data());
//	//	if constexpr (DEBUG) image_real<T, false>(p_n0_, p_n1_, temp1_r.data(), "6_border_conv.png");
//	//	//--------- Rectified solution & border --------
//	//	for (int32 n0 = 0; n0 < n0_; ++n0) {
//	//		for (int32 n1 = 0; n1 < n1_; ++n1) {
//	//			const int32 i1 = (n1 + off1_) + (n0 + off0_) * p_n1_;
//	//			const int32 i2 = (n1 + off1_ + h2) + (n0 + off0_ + w2) * p_n1_;
//	//			const int32 v = int32(std::abs((temp1_r[i2] * SCALE)));
//	//			temp2_r[i1] = temp2_r[i1] * v;
//	//		}
//	//	}
//	//	if constexpr (DEBUG) image_real<T, false>(p_n0_, p_n1_, temp2_r.data(), "7_border_Sol.png");
//
//	//};
//
//	std::vector<std::tuple<FIntVector, T, int32>> out;
//	//if (_perm | 0x1) {
//	//	perm(_ext[0], _ext[1], 0);
//	//	for (int32 n0 = 0; n0 < n0_; ++n0) {
//	//		for (int32 n1 = 0; n1 < n1_; ++n1) {
//	//			const int32 i1 = n1 + n0 * n1_;
//	//			const int32 i2 = (n1 + off1_) + (n0 + off0_) * p_n1_;
//	//			const T v = temp2_r[i2];
//	//			if (FMath::IsNearlyZero(v)) continue;
//	//			if (zh + map[i1] >= h_) continue;
//	//			out.push_back({ { n0, n1, map[i1]}, costf(n0, n1, map[i1], v), 0});
//	//		}
//	//	}
//	//}
//
//	//if (_perm | 0x2) {
//	//	perm(_ext[1], _ext[0], 1);
//	//	for (int32 n0 = 0; n0 < n0_; ++n0) {
//	//		for (int32 n1 = 0; n1 < n1_; ++n1) {
//	//			const int32 i1 = n1 + n0 * n1_;
//	//			const int32 i2 = (n1 + off1_) + (n0 + off0_) * p_n1_;
//	//			const T v = temp2_r[i2];
//	//			if (FMath::IsNearlyZero(v)) continue;
//	//			if (zh + map[i1] >= h_) continue;
//	//			out.push_back({ { n0, n1, map[i1] }, costf(n0, n1, map[i1], v), 1});
//	//		}
//	//	}
//	//}
//
//	//if (_perm | 0x8) {
//	//	perm(_ext[0], _ext[2], 2);
//	//	for (int32 n0 = 0; n0 < n0_; ++n0) {
//	//		for (int32 n1 = 0; n1 < n1_; ++n1) {
//	//			const int32 i1 = n1 + n0 * n1_;
//	//			const int32 i2 = (n1 + off1_) + (n0 + off0_) * p_n1_;
//	//			const T v = temp2_r[i2];
//	//			if (FMath::IsNearlyZero(v)) continue;
//	//			if (zh + map[i1] >= h_) continue;
//	//			out.push_back({ { n0, n1, map[i1] }, costf(n0, n1, map[i1], v), 2 });
//	//		}
//	//	}
//	//}
//
//	//if (_perm | 0x10) {
//	//	perm(_ext[2], _ext[0], 3);
//	//	for (int32 n0 = 0; n0 < n0_; ++n0) {
//	//		for (int32 n1 = 0; n1 < n1_; ++n1) {
//	//			const int32 i1 = n1 + n0 * n1_;
//	//			const int32 i2 = (n1 + off1_) + (n0 + off0_) * p_n1_;
//	//			const T v = temp2_r[i2];
//	//			if (FMath::IsNearlyZero(v)) continue;
//	//			if (zh + map[i1] >= h_) continue;
//	//			out.push_back({ { n0, n1, map[i1] }, costf(n0, n1, map[i1], v), 3 });
//	//		}
//	//	}
//	//}
//
//	//if (_perm | 0x20) {
//	//	perm(_ext[1], _ext[2], 4);
//	//	for (int32 n0 = 0; n0 < n0_; ++n0) {
//	//		for (int32 n1 = 0; n1 < n1_; ++n1) {
//	//			const int32 i1 = n1 + n0 * n1_;
//	//			const int32 i2 = (n1 + off1_) + (n0 + off0_) * p_n1_;
//	//			const T v = temp2_r[i2];
//	//			if (FMath::IsNearlyZero(v)) continue;
//	//			if (zh + map[i1] >= h_) continue;
//	//			out.push_back({ { n0, n1, map[i1] }, costf(n0, n1, map[i1], v), 4 });
//	//		}
//	//	}
//	//}
//
//	//if (_perm | 0x40) {
//	//	perm(_ext[2], _ext[1], 5);
//	//	for (int32 n0 = 0; n0 < n0_; ++n0) {
//	//		for (int32 n1 = 0; n1 < n1_; ++n1) {
//	//			const int32 i1 = n1 + n0 * n1_;
//	//			const int32 i2 = (n1 + off1_) + (n0 + off0_) * p_n1_;
//	//			const T v = temp2_r[i2];
//	//			if (FMath::IsNearlyZero(v)) continue;
//	//			if (zh + map[i1] >= h_) continue;
//	//			out.push_back({ { n0, n1, map[i1] }, costf(n0, n1, map[i1], v), 5 });
//	//		}
//	//	}
//	//}
//
//	return out;
//
//}//TurboPacker::Spectral::HeightMap::overlap
//
//template<class T, bool SIMD, bool DEBUG>
//void TurboPacker::Spectral::HeightMap<T, SIMD, DEBUG>::push(const FVector& _pos, const FVector& _ext) {
//
//	const int32 minY = int32(_pos[0] - _ext[0])-1;
//	const int32 maxY = int32(_pos[0] + _ext[0])+1;
//
//	const int32 minX = int32(_pos[1] - _ext[1])-1;
//	const int32 maxX = int32(_pos[1] + _ext[1])+1;
//
//	const int32 maxZ = int32(_pos[2] + _ext[2])+1;
//
//	//std::cout << "push: " << _pos << _ext << std::endl;
//
//	for (int32 n0 = minY; n0 < maxY; ++n0) {
//		for (int32 n1 = minX; n1 < maxX; ++n1) {
//			if (n0 < 0 || n0 >= config_->n0_ || n1 < 0 || n1 >= config_->n1_) {
//				//std::cout << "po: " << n0 << ", " << n1 << std::endl;
//				continue;
//			}
//			const int32 i = n1 + n0 * config_->n1_;
//			map[i] = std::max<T>(maxZ, map[i]);
//		}
//	}
//
//}//TurboPacker::Spectral::HeightMap::push
//
//template<class T, bool SIMD, bool DEBUG>
//void TurboPacker::Spectral::HeightMap<T, SIMD, DEBUG> ::print_size_in_bytes() const {
//	/*const int32 s1 = temp1_r.size() + temp2_r.size();
//	const int32 s2 = signal_c.size() + signal_log_c.size() + kernel_c.size() + sobel_x_c.size() + sobel_y_c.size() + border_c.size() + temp_c.size();
//	std::cout << ((s1 * sizeof(T) + s2 * sizeof(std::complex<T>)) /1000000) << "mb" << std::endl;*/
//}//TurboPacker::Spectral::HeightMap::print_size_in_bytes
//
//template<class T, bool SIMD, bool DEBUG>
//void TurboPacker::Spectral::HeightMap<T, SIMD, DEBUG>::save_heightmap(int32 _id) {
//	/*std::fill(temp1_r.begin(), temp1_r.end(), 0);
//	for (int32 n0 = 0; n0 < n0_; ++n0) {
//		for (int32 n1 = 0; n1 < n1_; ++n1) {
//			const int32 i1 = n1 + n0 * n1_;
//			const int32 i2 = (n1 + off1_) + (n0 + off0_) * p_n1_;
//			temp1_r[i2] = T(map[i1]);
//		}
//	}
//	for (int32 n0 = -1; n0 <= n0_; ++n0) {
//		const int32 i1 = (-1 + off1_) + (n0 + off0_) * p_n1_;
//		const int32 i2 = (n1_ + off1_) + (n0 + off0_) * p_n1_;
//		temp1_r[i1] = T(h_);
//		temp1_r[i2] = T(h_);
//	}
//	for (int32 n1 = -1; n1 <= n1_; ++n1) {
//		const int32 i1 = (n1 + off1_) + (-1 + off0_) * p_n1_;
//		const int32 i2 = (n1 + off1_) + (n0_ + off0_) * p_n1_;
//		temp1_r[i1] = T(h_);
//		temp1_r[i2] = T(h_);
//	}
//	std::stringstream ss;
//	ss << "aaa_" << _id << "_hm.png";
//	Debug::image_real<T, false>(p_n0_, p_n1_, temp1_r.data(), ss.str());*/
//}//TurboPacker::Spectral::HeightMap::save_heightmap
//
////----------------------------------------
////----------------------------------------
////----------------------------------------
//
//template<class C>
//void TurboPacker::Spectral::Debug::image_complex(
//	const int32 _n0,
//	const int32 _n1,
//	C _buffer,
//	const std::string& _filename,
//	double _max
//) {
//	{
//		const std::filesystem::path path = std::filesystem::path(TCHAR_TO_UTF8(*FPaths::ProjectDir())) / "overlap";
//		if (!std::filesystem::exists(path))
//			std::filesystem::create_directory(path);
//	}
//
//	std::vector<unsigned char> img;
//	for (int32 n0 = 0; n0 < _n0; ++n0) {
//		for (int32 n1 = 0; n1 < _n1; ++n1) {
//			const int32 i1 = n1 + n0 * _n1;
//
//			const std::complex<double> c(_buffer[i1].real(), _buffer[i1].imag());
//			const FVector col = Util::HSL_to_RGB_rad<double>(Util::c_to_HSL<double>(c, _max));
//
//			img.push_back(255 - (unsigned char)(col[0] * 255));
//			img.push_back(255 - (unsigned char)(col[1] * 255));
//			img.push_back(255 - (unsigned char)(col[2] * 255));
//			img.push_back(255);
//		}
//	}
//	const std::filesystem::path path = std::filesystem::path(TCHAR_TO_UTF8(*FPaths::ProjectDir())) / "overlap" / _filename;
//	lodepng::encode(path.string(), img.data(), _n1, _n0);
//
//}//TurboPacker::Spectral::Debug::image_complex
//
//template<class R, bool LogScaling>
//void TurboPacker::Spectral::Debug::image_real(
//	const int32 _n0,
//	const int32 _n1,
//	R* _buffer,
//	const std::string& _filename
//) {
//	{
//		const std::filesystem::path path = std::filesystem::path(TCHAR_TO_UTF8(*FPaths::ProjectDir())) / "overlap";
//		if (!std::filesystem::exists(path))
//			std::filesystem::create_directory(path);
//	}
//
//	double max = -std::numeric_limits<double>::infinity();
//	for (int32 i = 0; i < _n0 * _n1; ++i)
//		max = std::max<double>(max, std::abs(_buffer[i]));
//
//	std::vector<unsigned char> img;
//	for (int32 n0 = 0; n0 < _n0; ++n0) {
//		for (int32 n1 = 0; n1 < _n1; ++n1) {
//			const int32 i1 = n1 + n0 * _n1;
//			unsigned char c = 0;
//
//			if constexpr (LogScaling) {
//				c = std::abs(1. - std::abs(std::log(std::abs(1. + _buffer[i1]))) / std::log(max));
//			} else {
//				const double frac = 1. / max;
//				c = 255 * (std::abs(_buffer[i1]) * frac);
//			}
//
//			img.push_back(c);
//			img.push_back(c);
//			img.push_back(c);
//			img.push_back(255);
//		}
//	}
//	const std::filesystem::path path = std::filesystem::path(TCHAR_TO_UTF8(*FPaths::ProjectDir())) / "overlap/" / _filename;
//	lodepng::encode(path.string(), img.data(), _n1, _n0);
//}
