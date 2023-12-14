
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
		map.push(p.GetCenter(), p.GetExtent());
	}

	{
		const FBox p(FVector(20, 20, 0), FVector(48, 80, 3));
		map.push(p.GetCenter(), p.GetExtent());
	}

	{
		const FBox p(FVector(53, 20, 0), FVector(80, 80, 1));
		map.push(p.GetCenter(), p.GetExtent());
	}

	const auto start = std::chrono::high_resolution_clock::now();
	const auto res = map.overlap(FVector(11.));
	const std::chrono::duration<double> ee = std::chrono::high_resolution_clock::now() - start;
	std::cout << "Overlap: " << ee.count() << "s" << std::endl;
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

	using Rollcage = HeightMap<double, false, true>;
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
		}
		std::cout << r << "[" << std::exp(r) << "]" << std::endl;
	};

	{
		std::cout << "Kernel:" << std::endl;

		kernel[0] = 1;
		kernel[1] = 1;
		kernel[2] = 1;
		std::cout << kernel[0] << " | " << kernel[1] << " | " << kernel[2] << std::endl;

		kernel[3] = 1;
		kernel[4] = 1;
		kernel[5] = 1;
		std::cout << kernel[3] << " | " << kernel[4] << " | " << kernel[5] << std::endl;

		kernel[6] = 1;
		kernel[7] = 1;
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
