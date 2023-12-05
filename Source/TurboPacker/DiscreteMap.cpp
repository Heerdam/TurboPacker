
#include "DiscreteMap.h"

#include "lode/lodepng.h"
#include <unsupported/Eigen/FFT>
#include <math.h>

void ASpectralTester::Clear() {
	using namespace Util;

	UWorld* world = GetWorld();
	if (!world) return;

	FlushPersistentDebugLines(world);
}//ASpectralTester::Clear()

void ASpectralTester::TestInfGrid() {

	using namespace Util;

	UWorld* world = GetWorld();
	if (!world) return;

	FlushPersistentDebugLines(world);

	std::random_device rd;
	Rand rand(rd());

	Eigen::Vector3d size(10., 10., 10.);

	InfBucketGrid<bool, FVector> map(size);

	Dist d(0., 50.);
	for (int32 i = 0; i < 100000; ++i) {
		map.push(true, FVector(d(rand), d(rand), d(rand)));
	}

	for (int32 z = 0; z < 5; ++z) {
		for (int32 y = 0; y < 5; ++y) {
			for (int32 x = 0; x < 5; ++x) {

				DrawDebugBox(world, FVector(x * 10.f, y * 10.f, z * 10.f), FVector(4.9f), FColor::Blue, true);

			}
		}
	}

	auto b1 = map.radial_search(FVector(25.f), 5.f);

	std::cout << b1.size() << std::endl;

	for (const auto t : b1) {
		const auto& [tt, p] = *t;
		DrawDebugPoint(world, p, 1.f, FColor::Red, true);
	}

}//ASpectralTester::TestInfGrid

void ASpectralTester::TestHeightMap() {

	using namespace Util;

	UWorld* world = GetWorld();
	if (!world) return;

	FlushPersistentDebugLines(world);

	//const FBox bounds(FVector(0.f), FVector(100, 100, 1));
	//DrawDebugBox(world, bounds.GetCenter(), bounds.GetExtent(), FColor::Blue, true);

	//1200*800*1700
	using Rollcage = HeightMap<double>;
	//const auto start1 = std::chrono::high_resolution_clock::now();
	Rollcage map(100, 100, 1);
	//const std::chrono::duration<double> ee1 = std::chrono::high_resolution_clock::now() - start1;
	//std::cout << "Plan: " << ee1.count() << "s" << std::endl;

	//map.print_size_in_bytes();

	const FBox p(FVector(40, 40, 0), FVector(60, 60, 1));
	//DrawDebugBox(world, p.GetCenter() + FVector(5, 5, 0), p.GetExtent(), FColor::Blue, true);

	map.push(p.GetCenter(), p.GetExtent());

	/*
	for (int32 y = 0; y < map.y_; ++y) {
		for (int32 x = 0; x < map.x_; ++x) {
			const int32 i = map[{x, y}];

			if(i == 0)
				DrawDebugPoint(world, FVector(x, y, 0), 1.f, FColor::Blue, true);
			else
				DrawDebugPoint(world, FVector(x, y, 0), 1.f, FColor::Green, true);
		}
	}
	*/

	if constexpr (true) {
		const auto start = std::chrono::high_resolution_clock::now();
		const auto res = map.overlap(FVector(11.));
		const std::chrono::duration<double> ee = std::chrono::high_resolution_clock::now() - start;
		std::cout << "Overlap: " << ee.count() << "s" << std::endl;

		std::vector<unsigned char> img;
		for(int32 n0 = 0; n0 < map.n0_; ++n0){
			for (int32 n1 = 0; n1 < map.n1_; ++n1) {
				const int32 i = n1 + n0 * map.n1_;
				switch (res[0][i]) {
					case -1:
					{
						img.push_back(255);
						img.push_back(0);
						img.push_back(0);
						img.push_back(255);
					}
					break;
					case 0:
					{
						img.push_back(0);
						img.push_back(255);
						img.push_back(0);
						img.push_back(255);
					}
					break;
					case 1:
					{
						img.push_back(0);
						img.push_back(0);
						img.push_back(255);
						img.push_back(255);
					}
					break;
				}

				//DrawDebugPoint(world, FVector(x, y, 0), 1.f, FColor(0, 0, int32(double(res[1][i]) * frac)), true);
			}
			lodepng::encode("C:\\Users\\Heerdam\\Desktop\\Coding\\ue5\\TurboPacker\\test_r.png", img.data(), map.n1_, map.n0_);
		}





		//double max = 0;
		//for (int32 i = 0; i < map.p_n0_ * map.p_n1_; ++i) {
		//	//max = std::max<double>(max, std::abs(1. - std::abs(std::log(std::abs(1. + res[i])))));
		//	max = std::max<double>(max, std::abs(res[i]));
		//}

		//if constexpr (true){
		//	std::vector<unsigned char> img;
		//	for (int32 i = 0; i < map.p_n0_ * map.p_n1_; ++i) {
		//		
		//		//const std::complex<double> c(res[i], 0);
		//		//std::cout << c << std::endl;
		//		//const FVector hsl = c_to_HSL<double>(c, 10);
		//		//const FVector col = HSL_to_RGB_rad<double>(hsl);
		//		//const double mag =  std::abs(1. - std::abs(std::log(std::abs(1. + res[i]))) / max);
		//		const double mag = std::abs(res[i]) / max;
		//		//std::cout << mag << std::endl;

		//		img.push_back(255 - (unsigned char)(mag * 255.));
		//		img.push_back(255 - (unsigned char)(mag * 255.));
		//		img.push_back(255 - (unsigned char)(mag * 255.));


		//		//img.push_back(255 - int32(col[0] * 255.));
		//		//img.push_back(255 - int32(col[1] * 255.));
		//		//img.push_back(255 - int32(col[2] * 255.));
		//		img.push_back(255);
		//	}

		//	lodepng::encode("C:\\Users\\Heerdam\\Desktop\\Coding\\ue5\\TurboPacker\\test_r.png", img.data(), map.p_n1_, map.p_n0_);
		//}

		//const double frac = 1. / 1210. * 255.;

		

	}

}//ASpectralTester::TestHeightMap

void ASpectralTester::TestFFTW() {

	using namespace Util;

	const int32 n0 = 100;
	const int32 n1 = 100;

	int n1_cplx = n1 / 2 + 1;
	//printf("n1_cplx = %d\n", n1_cplx);

	double* in = fftw_alloc_real(n0 * n1);
	fftw_complex* out = fftw_alloc_complex(n0 * n1_cplx);
	fftw_complex* ref_out = fftw_alloc_complex(n0 * n1_cplx);

	fftw_plan p = fftw_plan_dft_r2c_2d(n0, n1, in, out, FFTW_ESTIMATE);

	//input
	for (int32 y = 0; y < n0; ++y) {
		for (int32 x = 0; x < n1; ++x) {

			const int32 i = x + y * n1;

			if ((x >= 40 && x <= 60) && (y >= 40 && y <= 60)) {
				in[i] = 1.;
			} else {
				in[i] = 0.;
			}

		}
	}

	// manually compute DFT for reference
	/*int idx_k, idx_j;
	double phi;
	for (int k0 = 0; k0 < n0; ++k0) {
		for (int k1 = 0; k1 < n1_cplx; ++k1) {
			idx_k = k0 * n1_cplx + k1;

			ref_out[idx_k] = 0.0;

			for (int j0 = 0; j0 < n0; ++j0) {
				for (int j1 = 0; j1 < n1; ++j1) {
					idx_j = j0 * n1 + j1;

					phi = -2.0 * PI * (k0 * j0 / ((double)n0)
						+ k1 * j1 / ((double)n1));

					ref_out[idx_k] += in[idx_j] * std::exp(I * phi);
				}
			}
		}
	}*/

	fftw_execute(p);

	{
		std::vector<unsigned char> img;
		for (int32 i = 0; i < n0 * n1; ++i) {
			img.push_back(int32(in[i] * 255.));
			img.push_back(int32(in[i] * 255.));
			img.push_back(int32(in[i] * 255.));
			img.push_back(255);
		}

		lodepng::encode("C:\\Users\\Heerdam\\Desktop\\Coding\\ue5\\TurboPacker\\test_r_ref.png", img.data(), n0, n1);
	}
	{
		std::vector<unsigned char> img;
		for (int32 i = 0; i < n0 * n1_cplx; ++i) {
			const std::complex<double> c(out[i][0], out[i][1]);
			//std::cout << c << std::endl;
			const FVector hsl = c_to_HSL<double>(c, 2.);
			const FVector col = HSL_to_RGB_rad<double>(hsl);
			//const double mag =  std::abs(1. - std::abs(std::log(std::abs(c))) / 10.);
			//std::cout << mag << std::endl;

			img.push_back(255 - int32(col[0] * 255.));
			img.push_back(255 - int32(col[1] * 255.));
			img.push_back(255 - int32(col[2] * 255.));
			img.push_back(255);
		}

		lodepng::encode("C:\\Users\\Heerdam\\Desktop\\Coding\\ue5\\TurboPacker\\test_c_ref.png", img.data(), n1_cplx, n0);
	}

	fftw_destroy_plan(p);
	fftw_free(in);
	fftw_free(out);
	fftw_free(ref_out);

}//ASpectralTester::TestFFTW