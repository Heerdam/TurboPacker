
#include "DiscreteMap.h"

#include "lode/lodepng.h"
#include <unsupported/Eigen/FFT>
#include <math.h>

FTransform TurboPacker::Detail::make_transform(
	EAxisPerm _perm,
	const FBox& _target,
	const FVector& _pivot_offset,
	const FVector& _relative_offset
) {
	const FVector& ctr = _target.GetCenter();

	switch (_perm) {
		case EAxisPerm::Z_XY_0:
		{
			FTransform tr;
			tr.SetRotation(FQuat::MakeFromEuler(FVector(0., 0., 0)));
			tr.SetLocation(ctr - _pivot_offset - _relative_offset);
			return tr;
		}
		case EAxisPerm::Z_XY_1:
		{
			FTransform tr;
			tr.SetRotation(FQuat::MakeFromEuler(FVector(0., 0., 90.)));
			tr.SetLocation(ctr + FVector(-_pivot_offset.Y, _pivot_offset.X, -_pivot_offset.Z) + FVector(_relative_offset.Y, -_relative_offset.X, -_relative_offset.Z));
			return tr;
		}
		case EAxisPerm::Z_XY_2:
		{
			FTransform tr;
			tr.SetRotation(FQuat::MakeFromEuler(FVector(0., 0., 2 * 90.)));
			tr.SetLocation(ctr + FVector(_pivot_offset.X, _pivot_offset.Y, -_pivot_offset.Z) + FVector(_relative_offset.X, _relative_offset.Y, -_relative_offset.Z));
			return tr;
		}
		case EAxisPerm::Z_XY_3:
		{
			FTransform tr;
			tr.SetRotation(FQuat::MakeFromEuler(FVector(0., 0., 3 * 90.)));
			tr.SetLocation(ctr + FVector(_pivot_offset.Y, _pivot_offset.X, -_pivot_offset.Z) - FVector(_relative_offset.Y, -_relative_offset.X, _relative_offset.Z));
			return tr;
		}
		//---------------------
		case EAxisPerm::Y_XZ_0:
		{
			FTransform tr;
			tr.SetRotation(FQuat::MakeFromEuler(FVector(90., 0., 0.)));
			tr.SetLocation(ctr + FVector(-_pivot_offset.X, -_pivot_offset.Z, _pivot_offset.Y) + FVector(-_relative_offset.X, -_relative_offset.Z, _relative_offset.Y));
			return tr;
		}
		case EAxisPerm::Y_XZ_1:
		{
			FTransform tr;
			tr.SetRotation(FQuat::MakeFromEuler(FVector(90., 0., 90.)));
			tr.SetLocation(ctr + FVector(_pivot_offset.Z, -_pivot_offset.X, _pivot_offset.Y) + FVector(_relative_offset.Z, -_relative_offset.X, _relative_offset.Y));
			return tr;
		}
		case EAxisPerm::Y_XZ_2:
		{
			FTransform tr;
			tr.SetRotation(FQuat::MakeFromEuler(FVector(90., 0., 2 * 90.)));
			tr.SetLocation(ctr + FVector(_pivot_offset.X, _pivot_offset.Z, _pivot_offset.Y) + FVector(_relative_offset.X, _relative_offset.Z, _relative_offset.Y));
			return tr;
		}
		case EAxisPerm::Y_XZ_3:
		{
			FTransform tr;
			tr.SetRotation(FQuat::MakeFromEuler(FVector(90., 0., 3 * 90.)));
			tr.SetLocation(ctr + FVector(-_pivot_offset.Z, _pivot_offset.X, _pivot_offset.Y) + FVector(-_relative_offset.Z, _relative_offset.X, _relative_offset.Y));
			return tr;
		}
		//---------------------
		case EAxisPerm::X_YZ_0:
		{
			FTransform tr;
			tr.SetRotation(FQuat::MakeFromEuler(FVector(90., 90., 0.)));
			tr.SetLocation(ctr + FVector(-_pivot_offset.Y, -_pivot_offset.Z, -_pivot_offset.X) + FVector(-_relative_offset.Y, -_relative_offset.Z, -_relative_offset.X));
			return tr;
		}
		case EAxisPerm::X_YZ_1:
		{
			FTransform tr;
			tr.SetRotation(FQuat::MakeFromEuler(FVector(2 * 90., 90., 0.)));
			tr.SetLocation(ctr + FVector(-_pivot_offset.Z, _pivot_offset.Y, -_pivot_offset.X) + FVector(-_relative_offset.Z, _relative_offset.Y, -_relative_offset.X));
			return tr;
		}
		case EAxisPerm::X_YZ_2:
		{
			FTransform tr;
			tr.SetRotation(FQuat::MakeFromEuler(FVector(3 * 90., 90., 0.)));
			tr.SetLocation(ctr + FVector(_pivot_offset.Y, _pivot_offset.Z, -_pivot_offset.X) + FVector(_relative_offset.Y, _relative_offset.Z, -_relative_offset.X));
			return tr;
		}
		case EAxisPerm::X_YZ_3:
		{
			FTransform tr;
			tr.SetRotation(FQuat::MakeFromEuler(FVector(0., 90., 0.)));
			tr.SetLocation(ctr + FVector(_pivot_offset.Z, -_pivot_offset.Y, -_pivot_offset.X) + FVector(_relative_offset.Z, -_relative_offset.Y, -_relative_offset.X));
			return tr;
		}
		//-------------------------
		//-------------------------
		case EAxisPerm::Z_n_XY_0:
		{
			FTransform tr;
			tr.SetRotation(FQuat::MakeFromEuler(FVector(0., 180., 0)));
			tr.SetLocation(ctr - FVector(-_pivot_offset.X, _pivot_offset.Y, -_pivot_offset.Z) + FVector(_relative_offset.X, -_relative_offset.Y, _relative_offset.Z));
			return tr;
		}
		case EAxisPerm::Z_n_XY_1:
		{
			FTransform tr;
			tr.SetRotation(FQuat::MakeFromEuler(FVector(0., 180., 90.)));
			tr.SetLocation(ctr + FVector(-_pivot_offset.Y, _pivot_offset.X, _pivot_offset.Z) + FVector(_relative_offset.Y, _relative_offset.X, _relative_offset.Z));
			return tr;
		}
		case EAxisPerm::Z_n_XY_2:
		{
			FTransform tr;
			tr.SetRotation(FQuat::MakeFromEuler(FVector(0., 180., 2 * 90.)));
			tr.SetLocation(ctr + FVector(-_pivot_offset.X, _pivot_offset.Y, _pivot_offset.Z) + FVector(-_relative_offset.X, _relative_offset.Y, _relative_offset.Z));
			return tr;
		}
		case EAxisPerm::Z_n_XY_3:
		{
			FTransform tr;
			tr.SetRotation(FQuat::MakeFromEuler(FVector(0., 180., 3 * 90.)));
			tr.SetLocation(ctr + FVector(-_pivot_offset.Y, -_pivot_offset.X, _pivot_offset.Z) + FVector(-_relative_offset.Y, -_relative_offset.X, _relative_offset.Z));
			return tr;
		}
		//-------------------------
		case EAxisPerm::Y_n_XZ_0:
		{
			FTransform tr;
			tr.SetRotation(FQuat::MakeFromEuler(FVector(90., 180., 0.)));
			tr.SetLocation(ctr + FVector(_pivot_offset.X, -_pivot_offset.Z, _pivot_offset.Y) + FVector(_relative_offset.X, -_relative_offset.Z, -_relative_offset.Y));
			return tr;
		}
		case EAxisPerm::Y_n_XZ_1:
		{
			FTransform tr;
			tr.SetRotation(FQuat::MakeFromEuler(FVector(90., 180., 90.)));
			tr.SetLocation(ctr + FVector(_pivot_offset.Z, _pivot_offset.X, _pivot_offset.Y) + FVector(_relative_offset.Z, _relative_offset.X, -_relative_offset.Y));
			return tr;
		}
		case EAxisPerm::Y_n_XZ_2:
		{
			FTransform tr;
			tr.SetRotation(FQuat::MakeFromEuler(FVector(90., 180., 2 * 90.)));
			tr.SetLocation(ctr + FVector(-_pivot_offset.X, _pivot_offset.Z, _pivot_offset.Y) + FVector(-_relative_offset.X, _relative_offset.Z, -_relative_offset.Y));
			return tr;
		}
		case EAxisPerm::Y_n_XZ_3:
		{
			FTransform tr;
			tr.SetRotation(FQuat::MakeFromEuler(FVector(90., 180., 3 * 90.)));
			tr.SetLocation(ctr + FVector(-_pivot_offset.Z, -_pivot_offset.X, _pivot_offset.Y) + FVector(-_relative_offset.Z, -_relative_offset.X, -_relative_offset.Y));
			return tr;
		}
		//-------------------------
		case EAxisPerm::X_n_YZ_0:
		{
			FTransform tr;
			tr.SetRotation(FQuat::MakeFromEuler(FVector(90., 180 + 90., 0.)));
			tr.SetLocation(ctr + FVector(-_pivot_offset.Y, -_pivot_offset.Z, _pivot_offset.X) + FVector(_relative_offset.Y, -_relative_offset.Z, _relative_offset.X));
			return tr;
		}
		case EAxisPerm::X_n_YZ_1:
		{
			FTransform tr;
			tr.SetRotation(FQuat::MakeFromEuler(FVector(90., -90., -90.)));
			tr.SetLocation(ctr + FVector(-_pivot_offset.Z, _pivot_offset.Y, _pivot_offset.X) + FVector(-_relative_offset.Z, -_relative_offset.Y, _relative_offset.X));
			return tr;
		}
		case EAxisPerm::X_n_YZ_2:
		{
			FTransform tr;
			tr.SetRotation(FQuat::MakeFromEuler(FVector(3 * 90., 180 + 90., 0.)));
			tr.SetLocation(ctr + FVector(_pivot_offset.Y, _pivot_offset.Z, _pivot_offset.X) + FVector(-_relative_offset.Y, _relative_offset.Z, _relative_offset.X));
			return tr;
		}
		case EAxisPerm::X_n_YZ_3:
		{
			FTransform tr;
			tr.SetRotation(FQuat::MakeFromEuler(FVector(-90., 270., 270.)));
			tr.SetLocation(ctr + FVector(_pivot_offset.Z, -_pivot_offset.Y, _pivot_offset.X) + FVector(_relative_offset.Z, _relative_offset.Y, _relative_offset.X));
			return tr;
		}
	}
	return FTransform();
}//TurboPacker::Detail::get_transform

//------------------------------------------

void ASpectralPacker::Pack() {
	using namespace TurboPacker;
	using namespace Spectral;

	UWorld* world = GetWorld();
	if (!world) return;

	Clear();

	HeightMap<double> map(Bounds.X, Bounds.Y, Bounds.Z);

	while (true) {

	}

}//ASpectralPacker::Pack

void ASpectralPacker::Clear() {
	UWorld* world = GetWorld();
	if (!world) return;

	FlushPersistentDebugLines(world);

	TArray<AActor*> tod;
	UGameplayStatics::GetAllActorsOfClass(world, APackerBox::StaticClass(), tod);
	for (AActor* a : tod)
		world->DestroyActor(a);
}//ASpectralPacker::Clear

//------------------------------------------

void ASpectralTester::Clear() {
	using namespace Util;

	UWorld* world = GetWorld();
	if (!world) return;

	FlushPersistentDebugLines(world);
}//ASpectralTester::Clear()

void ASpectralTester::TestHeightMap() {

	using namespace Util;
	using namespace TurboPacker;
	using namespace Spectral;

	UWorld* world = GetWorld();
	if (!world) return;

	Clear();

	//const FBox bounds(FVector(0.f), FVector(100, 100, 1));
	//DrawDebugBox(world, bounds.GetCenter(), bounds.GetExtent(), FColor::Blue, true);

	//1200*800*1700
	using Rollcage = HeightMap<double>;
	//const auto start1 = std::chrono::high_resolution_clock::now();
	Rollcage map(100, 100, 10);
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
		const auto res = map.overlap<true>(FVector(11.));
		const std::chrono::duration<double> ee = std::chrono::high_resolution_clock::now() - start;
		std::cout << "Overlap: " << ee.count() << "s" << std::endl;

		//std::vector<unsigned char> img;
		//for(int32 n0 = 0; n0 < map.n0_; ++n0){
		//	for (int32 n1 = 0; n1 < map.n1_; ++n1) {
		//		const int32 i = n1 + n0 * map.n1_;
		//		switch (res[0][i]) {
		//			case -1:
		//			{
		//				img.push_back(255);
		//				img.push_back(0);
		//				img.push_back(0);
		//				img.push_back(255);
		//			}
		//			break;
		//			case 0:
		//			{
		//				img.push_back(0);
		//				img.push_back(255);
		//				img.push_back(0);
		//				img.push_back(255);
		//			}
		//			break;
		//			case 1:
		//			{
		//				img.push_back(0);
		//				img.push_back(0);
		//				img.push_back(255);
		//				img.push_back(255);
		//			}
		//			break;
		//		}

		//		//DrawDebugPoint(world, FVector(x, y, 0), 1.f, FColor(0, 0, int32(double(res[1][i]) * frac)), true);
		//	}
		//	lodepng::encode("C:\\Users\\Heerdam\\Desktop\\Coding\\ue5\\TurboPacker\\test_r.png", img.data(), map.n1_, map.n0_);
		//}





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

void ASpectralTester::TestKernel() {

	std::array<int32, 9> test;
	std::array<int32, 9> kernel;

	const auto count = [&]() {
		int32 r = 0;
		for (int32 i = 0; i < 9; ++i)
			r += test[i] * kernel[i];
		std::cout << r << std::endl;
	};

	{
		std::cout << "Kernel:" << std::endl;

		/*
		-2 -2  1
		-2  0  1
		 1  1  2
		*/

		kernel[0] = -2;
		kernel[1] = -2;
		kernel[2] = 1;
		std::cout << kernel[0] << " | " << kernel[1] << " | " << kernel[2] << std::endl;

		kernel[3] = -2;
		kernel[4] = 0;
		kernel[5] = 1;
		std::cout << kernel[3] << " | " << kernel[4] << " | " << kernel[5] << std::endl;

		kernel[6] = 1;
		kernel[7] = 1;
		kernel[8] = 2;
		std::cout << kernel[6] << " | " << kernel[7] << " | " << kernel[8] << std::endl;

	}

	{
		std::cout << "Domain 0:" << std::endl;

		test[0] = 1;
		test[1] = 1;
		test[2] = 1;
		std::cout << test[0] << " | " << test[1] << " | " << test[2] << std::endl;

		test[3] = 1;
		test[4] = 1;
		test[5] = 1;
		std::cout << test[3] << " | " << test[4] << " | " << test[5] << std::endl;

		test[6] = 1;
		test[7] = 1;
		test[8] = 1;
		std::cout << test[6] << " | " << test[7] << " | " << test[8] << std::endl;
	}

	count();

	{
		std::cout << "Domain 1:" << std::endl;

		test[0] = 3;
		test[1] = 2;
		test[2] = 1;
		std::cout << test[0] << " | " << test[1] << " | " << test[2]  << std::endl;

		test[3] = 3;
		test[4] = 2;
		test[5] = 1;
		std::cout << test[3] << " | " << test[4] << " | " << test[5]  << std::endl;

		test[6] = 3;
		test[7] = 2;
		test[8] = 1;
		std::cout << test[6] << " | " << test[7] << " | " << test[8]  << std::endl;
	}

	count();

	{
		std::cout << "Domain 2:" << std::endl;

		test[0] = 1;
		test[1] = 2;
		test[2] = 3;
		std::cout << test[0] << " | " << test[1] << " | " << test[2] << std::endl;

		test[3] = 1;
		test[4] = 2;
		test[5] = 3;
		std::cout << test[3] << " | " << test[4] << " | " << test[5] << std::endl;

		test[6] = 1;
		test[7] = 2;
		test[8] = 3;
		std::cout << test[6] << " | " << test[7] << " | " << test[8] << std::endl;
	}

	count();

	{
		std::cout << "Domain 3:" << std::endl;

		test[0] = 3;
		test[1] = 3;
		test[2] = 3;
		std::cout << test[0] << " | " << test[1] << " | " << test[2] << std::endl;

		test[3] = 2;
		test[4] = 2;
		test[5] = 2;
		std::cout << test[3] << " | " << test[4] << " | " << test[5] << std::endl;

		test[6] = 1;
		test[7] = 1;
		test[8] = 1;
		std::cout << test[6] << " | " << test[7] << " | " << test[8] << std::endl;
	}

	count();

	{
		std::cout << "Domain 4:" << std::endl;

		test[0] = 1;
		test[1] = 1;
		test[2] = 1;
		std::cout << test[0] << " | " << test[1] << " | " << test[2] << std::endl;

		test[3] = 2;
		test[4] = 2;
		test[5] = 2;
		std::cout << test[3] << " | " << test[4] << " | " << test[5] << std::endl;

		test[6] = 3;
		test[7] = 3;
		test[8] = 3;
		std::cout << test[6] << " | " << test[7] << " | " << test[8] << std::endl;
	}

	count();

	{
		std::cout << "Domain 5:" << std::endl;

		test[0] = 1;
		test[1] = 1;
		test[2] = 1;
		std::cout << test[0] << " | " << test[1] << " | " << test[2] << std::endl;

		test[3] = 1;
		test[4] = 1;
		test[5] = 1;
		std::cout << test[3] << " | " << test[4] << " | " << test[5] << std::endl;

		test[6] = 2;
		test[7] = 2;
		test[8] = 2;
		std::cout << test[6] << " | " << test[7] << " | " << test[8] << std::endl;
	}

	count();

	{
		std::cout << "Domain 6:" << std::endl;

		test[0] = 2;
		test[1] = 1;
		test[2] = 1;
		std::cout << test[0] << " | " << test[1] << " | " << test[2] << std::endl;

		test[3] = 2;
		test[4] = 1;
		test[5] = 1;
		std::cout << test[3] << " | " << test[4] << " | " << test[5] << std::endl;

		test[6] = 2;
		test[7] = 1;
		test[8] = 1;
		std::cout << test[6] << " | " << test[7] << " | " << test[8] << std::endl;
	}

	count();

	{
		std::cout << "Domain 5:" << std::endl;

		test[0] = 2;
		test[1] = 2;
		test[2] = 2;
		std::cout << test[0] << " | " << test[1] << " | " << test[2] << std::endl;

		test[3] = 2;
		test[4] = 2;
		test[5] = 2;
		std::cout << test[3] << " | " << test[4] << " | " << test[5] << std::endl;

		test[6] = 1;
		test[7] = 1;
		test[8] = 1;
		std::cout << test[6] << " | " << test[7] << " | " << test[8] << std::endl;
	}

	count();

	{
		std::cout << "Domain 6:" << std::endl;

		test[0] = 1;
		test[1] = 2;
		test[2] = 2;
		std::cout << test[0] << " | " << test[1] << " | " << test[2] << std::endl;

		test[3] = 1;
		test[4] = 2;
		test[5] = 2;
		std::cout << test[3] << " | " << test[4] << " | " << test[5] << std::endl;

		test[6] = 1;
		test[7] = 2;
		test[8] = 2;
		std::cout << test[6] << " | " << test[7] << " | " << test[8] << std::endl;
	}

	count();

	{
		std::cout << "Domain 7:" << std::endl;

		test[0] = 2;
		test[1] = 2;
		test[2] = 2;
		std::cout << test[0] << " | " << test[1] << " | " << test[2] << std::endl;

		test[3] = 2;
		test[4] = 1;
		test[5] = 1;
		std::cout << test[3] << " | " << test[4] << " | " << test[5] << std::endl;

		test[6] = 2;
		test[7] = 1;
		test[8] = 1;
		std::cout << test[6] << " | " << test[7] << " | " << test[8] << std::endl;
	}

	count();

	{
		std::cout << "Domain 8:" << std::endl;

		test[0] = 1;
		test[1] = 1;
		test[2] = 1;
		std::cout << test[0] << " | " << test[1] << " | " << test[2] << std::endl;

		test[3] = 1;
		test[4] = 2;
		test[5] = 2;
		std::cout << test[3] << " | " << test[4] << " | " << test[5] << std::endl;

		test[6] = 1;
		test[7] = 2;
		test[8] = 2;
		std::cout << test[6] << " | " << test[7] << " | " << test[8] << std::endl;
	}

	count();

}//ASpectralTester::TestKernel