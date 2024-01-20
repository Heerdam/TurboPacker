
#include "PackTester.h"

#include "MQT.hpp"

FTransform Detail::make_transform(
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

//----------------------------------------

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

//----------------------------------------

void APackTester::Clear() {
	UWorld* world = GetWorld();
	if (!world) return;

	cache.Empty();
	idx = 0;
	
	FlushPersistentDebugLines(world);
	
	TArray<AActor*> tod;
	UGameplayStatics::GetAllActorsOfClass(world, APackerBox::StaticClass(), tod);
	for (AActor* a : tod)
		world->DestroyActor(a);
}

void APackTester::Pack() {
	UWorld* world = GetWorld();
	if (!world) return;

	Clear();

	using namespace MQT;

	std::vector<double> map;
	map.resize(Bounds.X * Bounds.Y);
	std::fill(map.begin(), map.end(), 0.);

	MedianQuadTree<double> tree(map, Bounds.X, Bounds.Y, Bounds.Z, BucketExtend);

	struct Result {
		double weight;
		int32 n0, n1, h;
		FVector ext;
		EAxisPerm perm;
		TSubclassOf<class APackerBox> box;
	};

	DrawDebugBox(world, 
		FVector(Bounds.X * 0.5, Bounds.Y * 0.5, Bounds.Z * 0.5), 
		FVector(Bounds.X * 0.5, Bounds.Y * 0.5, Bounds.Z * 0.5), FColor::Blue, true);

	int32 k = 0;
	const auto startt = std::chrono::high_resolution_clock::now();
	while (true) {
	//for(int32 k = 0; k < 5; ++k){

		std::vector<Result> res;

		std::cout << "----" << std::endl;

		for (const auto& b : Boxes) {

			const FBox aabb = b->GetDefaultObject<APackerBox>()->get_aabb();

			//Z_XY 
			const auto start = std::chrono::high_resolution_clock::now();
			for (int32 n0 = aabb.GetExtent().X; n0 < Bounds.X - aabb.GetExtent().X; ++n0) {
				for (int32 n1 = aabb.GetExtent().Y; n1 < Bounds.Y - aabb.GetExtent().Y; ++n1) {

					const int32 i = n1 + n0 * Bounds.Y;
					if (map[i] + aabb.GetSize().Z >= Bounds.Z) continue;

					const auto [l, m, h] = tree.check_overlap(
						Vec2{ int32_t(n0 - aabb.GetExtent().X), int32_t(n1 - aabb.GetExtent().Y) },
						Vec2{ int32_t(n0 + aabb.GetExtent().X), int32_t(n1 + aabb.GetExtent().Y) },
						map[i]);

					//const auto [l, m, h] = ::MQT::Detail::naive_tester<double>(
					//	map,
					//	Vec2{ int32_t(n0 - aabb.GetExtent().X), int32_t(n1 - aabb.GetExtent().Y) },
					//	Vec2{ int32_t(n0 + aabb.GetExtent().X), int32_t(n1 + aabb.GetExtent().Y) },
					//	Bounds.Y,
					//	map[i]);

					//if (l1 != l2 || m1 != m2 || h1 != h2) {
					//	std::cout << "[" << l1 << ", " << m1 << ", " << h1 << "]["
					//		<< "[" << l2 << ", " << m2 << ", " << h2 << "]" << std::endl;;
					//}

					//std::cout << l << ", " << m << ", " << h << std::endl;

					//if(k == 1)
						//std::cout << int32_t(n0 - aabb.GetExtent().X) << ", " << int32_t(n1 - aabb.GetExtent().Y) << "]["
							//<< int32_t(n0 + aabb.GetExtent().X) << ", " << int32_t(n1 + aabb.GetExtent().Y) << std::endl;

					if (h != 0 || l != 0) {
						//if (k == 1) {
						//	DrawDebugPoint(world, FVector(n0, n1, 0.), 2., FColor::Red, true);
							//DrawDebugBox(world,
								//FVector(double(n0), double(n1), 0.),
								//aabb.GetExtent(),
								//FColor::Red, true);
						//}
						continue;
					}

					res.emplace_back(
						
						std::pow(double(map[i]), 3) + std::pow(double(l), 2) + std::pow(n0 + n1, 2),
						n0, n1, map[i],
						aabb.GetExtent(),
						EAxisPerm::Z_XY_0, b
						
					);

					//if (k == 1) {
						//DrawDebugPoint(world, FVector(n0, n1, 0.), 2., FColor::Green, true);
						//DrawDebugBox(world,
						//	FVector(double(n0 - aabb.GetExtent().X), double(n1 - aabb.GetExtent().Y), 0.),
						//	FVector(double(n0 + aabb.GetExtent().X), double(n1 + aabb.GetExtent().Y), 1.),
						//	FColor::Green, true);
					//}
				}
			}

			const std::chrono::duration<double> ee = std::chrono::high_resolution_clock::now() - start;
			std::cout << ee.count() << std::endl;

		}

		//if (k == 1) return;

		if (res.empty()) break;

		std::sort(res.begin(), res.end(), [](const auto& _e1, const auto& _e2) {
			return _e1.weight < _e2.weight;
		});

		const auto& r = res.front();
		APackerBox* box = r.box->GetDefaultObject<APackerBox>();

		const FBox tar = FBox(
			FVector(r.n0 - r.ext.X, r.n1 - r.ext.Y, r.h), 
			FVector(r.n0 + r.ext.X, r.n1 + r.ext.Y, r.h + 2*r.ext.Z));

		DrawDebugBox(world, tar.GetCenter(), tar.GetExtent(), FColor::Red, true);

		FActorSpawnParameters params;
		params.bHideFromSceneOutliner = true;
		
		
		world->SpawnActor<AActor>(box->GetClass(), ::Detail::make_transform(
			r.perm,
			tar,
			box->get_aabb().GetCenter(),
			box->get_relative_location()
		), params);
		
		k++;
	
		for (int32 n0 = int32(tar.Min.X); n0 <= int32(tar.Max.X); ++n0) {
			for (int32 n1 = int32(tar.Min.Y); n1 <= int32(tar.Max.Y); ++n1) {
				const int32 i = n1 + n0 * Bounds.Y;

				map[i] = tar.Max.Z;

			}
		}

		tree.recompute();

		//std::stringstream ss;
		//ss << "map_" << k++ << ".png";
		//image_real<double, false>(Bounds.X, Bounds.Y, map.data(), ss.str());
	}

	const std::chrono::duration<double> eet = std::chrono::high_resolution_clock::now() - startt;
	std::cout << "Spawned " << k << " boxes in " << eet.count() << "s" << std::endl;
}

void APackTester::StepForward() {

}

void APackTester::StepBack() {

}

// Naive:
//  Spawned 162 boxes in 57.0335s - 0.352s /box
// MQT:
//  Spawned 162 boxes in 15.7559s - 0.09725s/box
//  Spawned 648 boxes in 160.478s - 0.247s/box