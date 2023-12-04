
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

	const FBox bounds(FVector(0.f), FVector(120, 80, 1700));
	//DrawDebugBox(world, bounds.GetCenter(), bounds.GetExtent(), FColor::Blue, true);

	//1200*800*1700
	using Rollcage = HeightMap<double>;
	const auto start1 = std::chrono::high_resolution_clock::now();
	Rollcage map(120, 80, 1700);
	const std::chrono::duration<double> ee1 = std::chrono::high_resolution_clock::now() - start1;
	std::cout << "Plan: " << ee1.count() << "s" << std::endl;

	map.print_size_in_bytes();

	const FBox p(FVector(40, 20, 0), FVector(80, 60, 10));
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
		const auto res = map.overlap(FVector(10.));
		const std::chrono::duration<double> ee = std::chrono::high_resolution_clock::now() - start;
		std::cout << "Overlap: " << ee.count() << "s" << std::endl;

		double max = -std::numeric_limits<double>::infinity();
		for (int32 i = 0; i < res.size(); ++i) {
			//std::cout << std::abs(res[1][i]) << std::endl;
			max = std::max(max, res[i].real());
		}

		std::cout << max << std::endl;

		std::vector<unsigned char> img;
		for (int32 y = 0; y < map.py_; ++y) {
			for (int32 x = 0; x < map.px_; ++x) {
				const int32 i = y + map.py_ * x;
				const auto c1 = c_to_HSL<double>(res[i], max);
				const auto col = FColor(int32(c1[0] / (2 * PI) * 255.));
				//const auto col = FColor(255 * int32(res[i].real() / max));
				img.push_back(col.R);
				img.push_back(col.G);
				img.push_back(col.B);
				img.push_back(255);
				//const auto c = HSL_to_RGB_rad<double>(c1);
				//std::cout << c1[0] << ", " << c1[1] << ", " << c1[2] << std::endl;
				//DrawDebugPoint(world, FVector(x, y, 0)*0.1, 2.5f, FColor(int32(c1[0] / (2*PI) * 255.)), true);
				//DrawDebugBox(world, FVector(x*10, y*10, 0), FVector(5., 5., 0.1), FColor(int32(c1[0] / (2 * PI) * 255.)), true);
			}
		}


		lodepng::encode("C:\\Users\\Heerdam\\Desktop\\Coding\\ue5\\TurboPacker\\test.png", img.data(), map.px_, map.py_);

		//const double frac = 1. / 1210. * 255.;

		//for(int32 x = 0; x < map.x_; ++x){
		//	for (int32 y = 0; y < map.y_; ++y) {
		//		const int32 i = map.kidx(x, y);
		//		switch (res[0][i]) {
		//			case -1:
		//			DrawDebugPoint(world, FVector(x, y, 0), 1.f, FColor::Blue , true);
		//			break;
		//			case 0:
		//			DrawDebugPoint(world, FVector(x, y, 0), 1.f, FColor::Green, true);
		//			break;
		//			case 1:
		//			DrawDebugPoint(world, FVector(x, y, 0), 1.f, FColor::Red, true);
		//			break;
		//		}

		//		//DrawDebugPoint(world, FVector(x, y, 0), 1.f, FColor(0, 0, int32(double(res[1][i]) * frac)), true);
		//	}
		//}

		

		//std::cout << max << std::endl;

		//for (int32 y = 0; y < map.y_; ++y) {
		//	for (int32 x = 0; x < map.x_; ++x) {
		//		const int32 i = map.kidx(x, y);
		//		//std::cout << res[1][i] << std::endl;

		//		//const int32 c = std::clamp(int32(double(std::abs(res[1][i])) / (double)(max / Itensity) * 255.), 0, 255);
		//		const int32 c = int32(((double(res[1][i]) - std::log(1)) * (double(max) - 1.)) / (std::log(max) - std::log(1)));
		//		DrawDebugPoint(world, FVector(x, y, 0), 1.5f, FColor(c, c, c), true);
		//		
		//		//DrawDebugPoint(world, FVector(x, y, 0), 1.f, FColor(255*res[1][i], 0, 0), true);
		//	}
		//}
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

		lodepng::encode("C:\\Users\\Heerdam\\Desktop\\Coding\\ue5\\TurboPacker\\test_r.png", img.data(), n0, n1);
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

		lodepng::encode("C:\\Users\\Heerdam\\Desktop\\Coding\\ue5\\TurboPacker\\test_c.png", img.data(), n1_cplx, n0);
	}

	fftw_destroy_plan(p);
	fftw_free(in);
	fftw_free(out);
	fftw_free(ref_out);

}//ASpectralTester::TestFFTW