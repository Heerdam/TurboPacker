
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

	HeightMap<double, false, true> map(Bounds.X, Bounds.Y, Bounds.Z);

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

	TArray<AActor*> tod;
	UGameplayStatics::GetAllActorsOfClass(world, APackerBox::StaticClass(), tod);
	for (AActor* a : tod)
		world->DestroyActor(a);
}//ASpectralTester::Clear()

//------------------------------------------

APackerBox::APackerBox() {
	PrimaryActorTick.bCanEverTick = false;
	RootComponent = CreateDefaultSubobject<USceneComponent>(TEXT("Root"));
	mesh = CreateDefaultSubobject<UStaticMeshComponent>(TEXT("Mesh"));
	mesh->SetupAttachment(RootComponent);
}//APackerBox::APackerBox

FBox APackerBox::get_aabb() const {
	const FVector s = mesh->GetRelativeScale3D();
	const FBox bb = mesh->GetStaticMesh()->GetBounds().GetBox();
	const FVector se = FVector(bb.GetExtent().X * s.X, bb.GetExtent().Y * s.Y, bb.GetExtent().Z * s.Z);

	const FBox out = FBox(bb.GetCenter() * s - se, bb.GetCenter() * s + se);

	const FRotator rot = mesh->GetRelativeRotation();
	const FVector re = rot.RotateVector(out.GetExtent());
	const FVector c = out.GetCenter();

	const FVector min = {
		std::min(c.X - re.X, c.X + re.X),
		std::min(c.Y - re.Y, c.Y + re.Y),
		std::min(c.Z - re.Z, c.Z + re.Z)
	};

	const FVector max = {
		std::max(c.X - re.X, c.X + re.X),
		std::max(c.Y - re.Y, c.Y + re.Y),
		std::max(c.Z - re.Z, c.Z + re.Z)
	};

	return FBox(min, max);
}//APackerBox::get_aabb

//------------------------------------------

void ASpectralTester::TestHeightMap() {
	using namespace Util;
	using namespace TurboPacker;
	using namespace Spectral;

	using Rollcage = HeightMap<double, false, true>;
	Rollcage map(100, 100, 10);

	map.print_size_in_bytes();

	{
		const FBox p(FVector(20, 20, 0), FVector(53, 80, 2));
		//map.push(p.GetCenter(), p.GetExtent());
	}

	{
		const FBox p(FVector(20, 20, 0), FVector(48, 80, 3));
		//map.push(p.GetCenter(), p.GetExtent());
	}

	{
		const FBox p(FVector(53, 20, 0), FVector(80, 80, 1));
		//map.push(p.GetCenter(), p.GetExtent());
	}

	const auto start = std::chrono::high_resolution_clock::now();
	const auto res = map.overlap(FVector(11.));
	const std::chrono::duration<double> ee = std::chrono::high_resolution_clock::now() - start;
	std::cout << "Overlap: " << ee.count() << "s" << std::endl;
}

void ASpectralTester::BenchConv() {
	using namespace Util;
	using namespace TurboPacker;
	using namespace Spectral;
	using namespace Spectral::Detail;

	UWorld* world = GetWorld();
	if (!world) return;

	{
		using Rollcage = HeightMap<double, false, false>;
		Rollcage map(1200, 1200, 1700);

		{
			const auto start = std::chrono::high_resolution_clock::now();
			map.overlap(FVector(100., 100., 100.));
			const std::chrono::duration<double> ee = std::chrono::high_resolution_clock::now() - start;
			std::cout << "Spectral 100: " << ee.count() << "s" << std::endl;
		}
		{
			const auto start = std::chrono::high_resolution_clock::now();
			map.overlap(FVector(700., 700, 100.));
			const std::chrono::duration<double> ee = std::chrono::high_resolution_clock::now() - start;
			std::cout << "Spectral 700: " << ee.count() << "s" << std::endl;
		}
	}

	{
		FFTWVector<double> map;
		map.resize(1200 * 1200);

		const auto idx = [](int32 _n0, int32 _n1) {
			return _n1 + _n0 * 1200;
		};

		int32 l = 0;
		int32 h = 0;

		const auto start = std::chrono::high_resolution_clock::now();
		for (int32 n0 = 49; n0 < 1200 - 49; ++n0) {
			for (int32 n1 = 49; n1 < 1200 - 49; ++n1) {

				const int32 i = idx(n0, n1);
				const int32 e = map[i];

				for (int m0 = -49; m0 < 50; ++m0) {
					for (int m1 = -49; m1 < 50; ++m1) {
						const int32 ii = idx(n0 + m0, n1 + m1);
						const int32 v = map[ii];
						if (v <= e) l++;
						else h++;
					}
				}
				
			}
		}
		std::cout << l << " " << h << std::endl;
		const std::chrono::duration<double> ee = std::chrono::high_resolution_clock::now() - start;
		std::cout << "naive 100: " << ee.count() << "s" << std::endl;

	}

	{
		FFTWVector<double> map;
		map.resize(1200 * 1200);

		const auto idx = [](int32 _n0, int32 _n1) {
			return _n1 + _n0 * 1200;
		};

		int32 l = 0;
		int32 h = 0;

		const auto start = std::chrono::high_resolution_clock::now();
		for (int32 n0 = 349; n0 < 1200 - 349; ++n0) {
			for (int32 n1 = 349; n1 < 1200 - 349; ++n1) {

				const int32 i = idx(n0, n1);
				const int32 e = map[i];

				for (int m0 = -349; m0 < 350; ++m0) {
					for (int m1 = -349; m1 < 350; ++m1) {
						const int32 ii = idx(n0 + m0, n1 + m1);
						const int32 v = map[ii];
						if (v <= e) l++;
						else h++;
					}
				}

			}
		}
		std::cout << l << " " << h << std::endl;
		const std::chrono::duration<double> ee = std::chrono::high_resolution_clock::now() - start;
		std::cout << "naive 700: " << ee.count() << "s" << std::endl;

	}
}

void ASpectralTester::TestPacker(){

	using namespace Util;
	using namespace TurboPacker;
	using namespace Spectral;

	UWorld* world = GetWorld();
	if (!world) return;

	if (Boxes.IsEmpty()) return;

	Clear();

	std::random_device rd;
	Rand r(rd());

	using Rollcage = HeightMap<double, false, false>;
	Rollcage map(Bounds.X, Bounds.Y, Bounds.Z);
	//Rollcage map(1200, 800, 1700);
	map.print_size_in_bytes();



	const FBox bounds(FVector(0.f), FVector(Bounds.X, Bounds.Y, Bounds.Z));
	DrawDebugBox(world, bounds.GetCenter(), bounds.GetExtent(), FColor::Blue, true);

	FActorSpawnParameters params;
	params.bHideFromSceneOutliner = true;
	params.SpawnCollisionHandlingOverride = ESpawnActorCollisionHandlingMethod::AlwaysSpawn;

	std::cout << "Packing started..." << std::endl;

	int32 count = 0;
	const auto start = std::chrono::high_resolution_clock::now();
	while (true) {
		//std::cout << "############### " << count << " ###############" << std::endl;
		//map.save_heightmap(count);
		//if (count == 10) break;
		
		const auto next = Boxes[Dist(0, Boxes.Num() - 1)(r)];

		const auto aabb = next->GetDefaultObject<APackerBox>()->get_aabb();
		//std::cout << "Next box: " << aabb << std::endl;
		const auto res = map.overlap(aabb.GetSize());
		//std::cout << "res: " << res.size() << std::endl;
		if (res.empty()) break;

		

		

		const auto rit = std::min_element(res.begin(), res.end(), [](const auto& _e1, const auto _e2) {
			const auto& [pos1, d1, p1] = _e1;
			const auto& [pos2, d2, p2] = _e2;
			return d2 > d1;
		});

		/*if (count == 2) {
			std::cout << "--------------" << std::endl;
			for (const auto& [pos, d, p] : res) {
				std::cout << "pos: " << pos << std::endl;
				std::cout << "p: " << p << std::endl;
				std::cout << "cost: " << d << std::endl;
			}
			std::cout << "--------------" << std::endl;
		}*/

		const auto& [pos, d, p] = *rit;
		//std::cout << "pos: " << pos << std::endl;
		//std::cout << "p: " << p << std::endl;
		//std::cout << "cost: " << d << std::endl;

		const FVector delta(pos.X, pos.Y, pos.Z);
		const FVector ext = aabb.GetExtent();

		switch (p) {
			case 0:
			{
				const FBox spos = FBox(FVector(-ext.X, -ext.Y, -ext.Z), FVector(ext.X, ext.Y, ext.Z)).ShiftBy(delta);
				map.push(spos.GetCenter() + FVector(0, 0, spos.GetExtent().Z), spos.GetExtent());
				//std::cout << "spos: "  << spos << std::endl;
				const auto trf = ::TurboPacker::Detail::make_transform(
					EAxisPerm::Z_XY_0,
					spos,
					aabb.GetCenter(),
					next->GetDefaultObject<APackerBox>()->get_relative_location()
				);
				world->SpawnActor<APackerBox>(next, trf, params);
			}
			break;
			case 1:
			{
				const FBox spos = FBox(FVector(-ext.X, -ext.Y, -ext.Z), FVector(ext.X, ext.Y, ext.Z)).ShiftBy(delta);
				map.push(spos.GetCenter() + FVector(0, 0, spos.GetExtent().Z), spos.GetExtent());
				//std::cout << spos << std::endl;
				const auto trf = ::TurboPacker::Detail::make_transform(
					EAxisPerm::Z_XY_1,
					spos,
					aabb.GetCenter(),
					next->GetDefaultObject<APackerBox>()->get_relative_location()
				);
				world->SpawnActor<APackerBox>(next, trf, params);
			}
			break;
			case 2:
			{
				const FBox spos = FBox(FVector(-ext.X, -ext.Y, -ext.Z), FVector(ext.X, ext.Y, ext.Z)).ShiftBy(delta);
				map.push(spos.GetCenter() + FVector(0, 0, spos.GetExtent().Z), spos.GetExtent());
				//std::cout << spos << std::endl;
				const auto trf = ::TurboPacker::Detail::make_transform(
					EAxisPerm::Y_XZ_0,
					spos,
					aabb.GetCenter(),
					next->GetDefaultObject<APackerBox>()->get_relative_location()
				);
				world->SpawnActor<APackerBox>(next, trf, params);
			}
			break;
			case 3:
			{
				const FBox spos = FBox(FVector(-ext.X, -ext.Y, -ext.Z), FVector(ext.X, ext.Y, ext.Z)).ShiftBy(delta);
				map.push(spos.GetCenter() + FVector(0, 0, spos.GetExtent().Z), spos.GetExtent());
				//std::cout << spos << std::endl;
				const auto trf = ::TurboPacker::Detail::make_transform(
					EAxisPerm::Y_XZ_1,
					spos,
					aabb.GetCenter(),
					next->GetDefaultObject<APackerBox>()->get_relative_location()
				);
				world->SpawnActor<APackerBox>(next, trf, params);
			}
			break;
			case 4:
			{
				const FBox spos = FBox(FVector(-ext.X, -ext.Y, -ext.Z), FVector(ext.X, ext.Y, ext.Z)).ShiftBy(delta);
				map.push(spos.GetCenter() + FVector(0, 0, spos.GetExtent().Z), spos.GetExtent());
				//std::cout << spos << std::endl;
				const auto trf = ::TurboPacker::Detail::make_transform(
					EAxisPerm::X_YZ_0,
					spos,
					aabb.GetCenter(),
					next->GetDefaultObject<APackerBox>()->get_relative_location()
				);
				world->SpawnActor<APackerBox>(next, trf, params);
			}
			break;
			case 5:
			{
				const FBox spos = FBox(FVector(-ext.X, -ext.Y, -ext.Z), FVector(ext.X, ext.Y, ext.Z)).ShiftBy(delta);
				map.push(spos.GetCenter() + FVector(0, 0, spos.GetExtent().Z), spos.GetExtent());
				//std::cout << spos << std::endl;
				const auto trf = ::TurboPacker::Detail::make_transform(
					EAxisPerm::X_YZ_1,
					spos,
					aabb.GetCenter(),
					next->GetDefaultObject<APackerBox>()->get_relative_location()
				);
				world->SpawnActor<APackerBox>(next, trf, params);
			}
			break;
		}

		count++;
		
	}
	const std::chrono::duration<double> ee = std::chrono::high_resolution_clock::now() - start;
	std::cout << "Packing: " << ee.count() << "s" << " [" << count << "]" << std::endl;

}//ASpectralTester::TestHeightMap

void ASpectralTester::TestKernel() {

	std::array<int32, 9> test;
	std::array<int32, 9> kernel;

	const auto count = [&]() {
		double r = 0;
		double r1 = 0;
		for (int32 i = 0; i < 9; ++i) {
			r += std::log(double(test[i] + 1)) * (double)kernel[i];
			r1 += test[i] * kernel[i];
		}
		std::cout << "lin: " << r1 << std::endl;
		std::cout << r << "[" << std::exp(r) << "]" << std::endl;
	};

	{
		std::cout << "Kernel:" << std::endl;

		kernel[0] = -1;
		kernel[1] = 0;
		kernel[2] = 1;
		std::cout << kernel[0] << " | " << kernel[1] << " | " << kernel[2] << std::endl;

		kernel[3] = -2;
		kernel[4] = 0;
		kernel[5] = 2;
		std::cout << kernel[3] << " | " << kernel[4] << " | " << kernel[5] << std::endl;

		kernel[6] = -1;
		kernel[7] = 0;
		kernel[8] = 1;
		std::cout << kernel[6] << " | " << kernel[7] << " | " << kernel[8] << std::endl;

	}

	{
		std::cout << "Domain 0:" << std::endl;

		test[0] = 2;
		test[1] = 2;
		test[2] = 2;
		std::cout << test[0] << " | " << test[1] << " | " << test[2] << std::endl;

		test[3] = 2;
		test[4] = 2;
		test[5] = 2;
		std::cout << test[3] << " | " << test[4] << " | " << test[5] << std::endl;

		test[6] = 2;
		test[7] = 2;
		test[8] = 2;
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

void ASpectralTester::QuadTreeTester() {

	using namespace Util;
	using namespace TurboPacker;
	using namespace Spectral;

	UWorld* world = GetWorld();
	if (!world) return;

	Clear();

	using Con = Spectral::Impl::Config<double, false, false>;

	constexpr int32 map0 = 30;
	constexpr int32 map1 = 30;
	constexpr int32 maph = 100;

	Con config (map0, map1, maph);
	config.world = world;

	 

	FFTWVector<Con::T> map;
	map.resize(map0 * map0);
	std::fill(map.begin(), map.end(), 0.);

	//-----------------------

	for (int32 n0 = 0; n0 < 10; ++n0) {
		for (int32 n1 = 0; n1 < 10; ++n1) {
			const int32 i = n1 + n0 * map1;
			map[i] = 1.;
		}
	}

	for (int32 n0 = 10; n0 < 20; ++n0) {
		for (int32 n1 = 10; n1 < 20; ++n1) {
			const int32 i = n1 + n0 * map1;
			map[i] = 2.;
		}
	}

	for (int32 n0 = 20; n0 < 30; ++n0) {
		for (int32 n1 = 20; n1 < 30; ++n1) {
			const int32 i = n1 + n0 * map1;
			map[i] = 3.;
		}
	}

	//-----------------------

	MedianQuadTree<Con::T, Con::SIMD, Con::DEBUG> tree (config, map, 10);
	tree.recompute();

	const auto[l, h] = tree.check_height(FIntVector(15, 15, 2), FIntVector2(10, 10));

	std::cout << "l: " << l << std::endl;
	std::cout << "h: " << h << std::endl;

}//ASpectralTester::QuadTreeTester

void ASpectralTester::QuadTreeTester2() {

	using namespace Util;
	using namespace TurboPacker;
	using namespace Spectral;

	UWorld* world = GetWorld();
	if (!world) return;

	Clear();

	//---------------------

	//Rollcage map(1200, 800, 1700);
	using Con = Spectral::Impl::Config<double, false, false>;

	constexpr int32 map0 = 800;
	constexpr int32 map1 = 1200;
	constexpr int32 maph = 1700;

	FFTWVector<Con::T> map;
	map.resize(map0 * map1);

	//---------------------

	const auto naive = [&](
		const int32 _p0,
		const int32 _p1,
		const int32 _ext,
		const int32 _h
		) {
			const int32 min_n0 = _p0 - _ext;
			const int32 max_n0 = _p0 + _ext;

			const int32 min_n1 = _p1 - _ext;
			const int32 max_n1 = _p1 + _ext;

			int32 h = 0;
			int32 l = 0;
			for (int32 n0 = min_n0; n0 < max_n0; ++n0) {
				for (int32 n1 = min_n1; n1 < max_n1; ++n1) {
					const int32 i = n1 + n0 * map1;
					if (map[i] >= _h) h++;
					else l++;
				}
			}

			return std::make_pair(l, h);
		};

	//---------------------

	Con config(map0, map1, maph);
	config.world = world;

	MedianQuadTree<Con::T, Con::SIMD, Con::DEBUG> tree(config, map, 10);

	for (int32 k = 0; k < 10; ++k) {
		std::cout << "--------------------------" << std::endl;
		std::cout << "complexity: " << k << std::endl;

		for (int32 n0 = k * 75; n0 < 75 + k * 75; ++n0) {
			for (int32 n1 = k * 75; n1 < 75 + k * 75; ++n1) {
				const int32 i = n1 + n0 * map1;
				map[i] = k * 8;
			}
		}

		tree.recompute();

		//---------------------
		const auto [l1, h1] = naive(map0 / 2, map1 / 2, 25 + k * 20, k * 5);
		const auto [l2, h2] = tree.check_height(FIntVector(map0 / 2, map1 / 2, k * 5), FIntVector2(25 + k * 20));
		
		if (l1 == l2 && h1 == h2) {
			std::cout << "correct" << std::endl;
		} else {
			std::cout << "false. exp: " << l1 << ", " << h1 << " - res: " << l2 << ", " << h2 << std::endl;
		}
	}

}//ASpectralTester::QuadTreeBench

void ASpectralTester::QuadTreeBench() {

	using namespace Util;
	using namespace TurboPacker;
	using namespace Spectral;

	UWorld* world = GetWorld();
	if (!world) return;

	Clear();

	//---------------------

	//Rollcage map(1200, 800, 1700);
	using Con = Spectral::Impl::Config<double, false, false>;

	constexpr int32 map0 = 800;
	constexpr int32 map1 = 1200;
	constexpr int32 maph = 1700;

	FFTWVector<Con::T> map;
	map.resize(map0 * map1);

	//---------------------

	const auto naive = [&](
		const int32 _p0,
		const int32 _p1,
		const int32 _ext,
		const int32 _h
		) {
			const int32 min_n0 = _p0 - _ext;
			const int32 max_n0 = _p0 + _ext;

			const int32 min_n1 = _p1 - _ext;
			const int32 max_n1 = _p1 + _ext;

			int32 h = 0;
			int32 l = 0;
			for (int32 n0 = min_n0; n0 < max_n0; ++n0) {
				for (int32 n1 = min_n1; n1 < max_n1; ++n1) {
					const int32 i = n1 + n0 * map1;
					if (map[i] >= _h) h++;
					else l++;
				}
			}

			return std::make_pair(l, h);
		};

	//---------------------

	Con config(map0, map1, maph);
	config.world = world;

	MedianQuadTree<Con::T, Con::SIMD, Con::DEBUG> tree(config, map, 10);

	for (int32 k = 0; k < 10; ++k) {
		std::cout << "--------------------------" << std::endl;
		std::cout << "complexity: " << k << std::endl;

		for (int32 n0 = k*75; n0 < 75 + k * 75; ++n0) {
			for (int32 n1 = k * 75; n1 < 75 + k * 75; ++n1) {
				const int32 i = n1 + n0 * map1;
				map[i] = k * 15;
			}
		}

		tree.recompute();

		//---------------------
		int32 x = 0;
		std::cout << "Naive" << std::endl;
		for (int32 i = 0; i < 10; ++i) {

			double t = 0;

			for (int32 j = 0; j < 15; ++j) {
				const auto start = std::chrono::high_resolution_clock::now();

				const auto [l, h] = naive(map0 / 2, map1 / 2, 25 + i * 20, i * 20);
				x += l + h;
				const std::chrono::duration<double> ee = std::chrono::high_resolution_clock::now() - start;
				t += ee.count();

				//std::cout << l << ", " << h << std::endl;
			}
			std::cout << i << ": " << t / 15. << "s" << std::endl;

		}

		//---------------------

		std::cout << "Quad" << std::endl;
		for (int32 i = 0; i < 10; ++i) {

			double t = 0;

			for (int32 j = 0; j < 15; ++j) {
				const auto start = std::chrono::high_resolution_clock::now();

				const auto [l, h] = tree.check_height(FIntVector(map0 / 2, map1 / 2, i * 20), FIntVector2(25 + i * 20));

				const std::chrono::duration<double> ee = std::chrono::high_resolution_clock::now() - start;
				t += ee.count();
				x += l + h;
				//std::cout << l << ", " << h << std::endl;
			}
			std::cout << i << ": " << t / 15. << "s" << std::endl;
		}

		std::cout << x << std::endl;
	}

}//ASpectralTester::QuadTreeBench

void ASpectralTester::BoxStackBench() {

	using namespace Util;
	using namespace TurboPacker;
	using namespace Spectral;

	UWorld* world = GetWorld();
	if (!world) return;

	Clear();

	//---------------------

	//Rollcage map(1200, 800, 1700);
	using Con = Spectral::Impl::Config<double, false, false>;

	constexpr int32 map0 = 800;
	constexpr int32 map1 = 1200;
	constexpr int32 maph = 1700;

	BoxStack<double> stack;

	std::random_device rd;
	Rand rand(rd());

	const auto rand_box = [&]() {

		const double s0 = DistD(100, 700)(rd);
		const double s1 = DistD(200, 1100)(rd);
		const double h = DistD(50, 1600)(rd);

		const double p0 = DistD(0, map0 - s0)(rd);
		const double p1 = DistD(0, map1 - s1)(rd);
		const double ph = DistD(0, maph - h)(rd);

		return FBox(FVector(p0, p1, ph), FVector(p0 + s0, p1 + s1, ph + h));
	};

	constexpr int32 sc = 25000;

	std::vector<FBox> test;
	test.reserve(sc);
	for (int32 i = 0; i < sc; ++i) {
		stack.push(rand_box());
		test.push_back(rand_box());
	}

	{
		int32 c = 0;
		const auto start = std::chrono::high_resolution_clock::now();

		for (int32 i = 0; i < sc; ++i) {
			if (stack.overlap(FBox(FVector(-50.), FVector(-10.))))
				c++;
		}

		const std::chrono::duration<double> ee = std::chrono::high_resolution_clock::now() - start;
		std::cout << ee.count() << std::endl;
		std::cout << c << std::endl;
	}

	{
		int32 c = 0;
		const auto start = std::chrono::high_resolution_clock::now();

		for (int32 i = 0; i < sc; ++i) {
			if (stack.overlap_naive(FBox(FVector(-50.), FVector(-10.))))
				c++;
		}

		const std::chrono::duration<double> ee = std::chrono::high_resolution_clock::now() - start;
		std::cout << ee.count() << std::endl;
		std::cout << c << std::endl;
	}

}//ASpectralTester::BoxStackBench