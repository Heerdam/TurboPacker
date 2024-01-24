
#include "PackTester.h"

#include "MQT.hpp"
#include "MQT2.hpp"

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

	using namespace MQT2;

	std::vector<double> map;
	map.resize(Bounds * Bounds);
	std::fill(map.begin(), map.end(), 0.);

	std::cout << Bounds << std::endl;

	using Tree = MedianQuadTree<double, 15>;

	Tree tree(map, Bounds);

	struct Result { 
		double weight;
		int32 n0, n1, h;
		FVector ext;
		EAxisPerm perm;
		TSubclassOf<class APackerBox> box;
	};

	DrawDebugBox(world, 
		FVector(Bounds * 0.5, Bounds * 0.5, Height * 0.5),
		FVector(Bounds * 0.5, Bounds * 0.5, Height * 0.5), FColor::Blue, true);

	int32 kkk = 0;
	int32 k = 0;
	const auto overlap = [&](const int32 _ext0, const int32 _ext1, const int32 _h)->std::optional<Result> {

		for (int32 n0 = _ext0; n0 < Bounds - _ext0; ++n0) {
			for (int32 n1 = _ext1; n1 < Bounds - _ext1; ++n1) {

				const int32 i = n1 + n0 * Bounds;
				//if (i < 0 || i >= map.size())
					//std::cout << "oob overlap" << std::endl;
				if (map[i] + 2 * _h >= Height) continue;

				const auto [l, m, h] = tree.check_overlap(
					Vec2{ int32_t(n0 - _ext0), int32_t(n1 - _ext1) },
					Vec2{ int32_t(n0 + _ext0), int32_t(n1 + _ext1) },
					map[i]);
				kkk++;

				
				/*const auto [l1, m1, h1] = ::MQT3::Detail::naive_tester<double>(
					map,
					Vec2{ int32_t(n0 - _ext0), int32_t(n1 - _ext1) },
					Vec2{ int32_t(n0 + _ext0), int32_t(n1 + _ext1) },
					Bounds,
					map[i]);
					
				

				if (l != l1 || m != m1 || h != h1) {
					std::cout << "-----------" << std::endl;
					std::cout << k << std::endl;
					std::cout << n0 << ", " << n1 << std::endl;
					std::cout << _ext0 << ", " << _ext1 << std::endl;
					std::cout << l << ", " << m << ", " << h << " | " << l1 << ", " << m1 << ", " << h1 << std::endl;
					DrawDebugBox(world, FVector(n0, n1, map[i] + _h), FVector(_ext0, _ext1, _h), FColor::Red, true);
				}*/

				if (h != 0 || l != 0) continue;

				Result out;
				out.n0 = n0;
				out.n1 = n1;
				out.h = map[i];
				out.weight = std::pow(double(map[i]), 3) + std::pow(double(l), 2) + std::pow(n0 + n1, 2);

				return { out };
			}
		}
		return {};
	};

	const int32_t bc = (Bounds / Tree::BUCKET_SIZE);
	std::vector<bool> mm;
	mm.resize(bc* bc);
	std::fill(mm.begin(), mm.end(), true);

	std::vector<std::pair<TSubclassOf<APackerBox>, FTransform>> toSpawn;
	 
	double vol = 0.;
	double s1 = 0., s2 = 0., s3 = 0., s4 = 0., s5 = 0.;
	const auto startt = std::chrono::high_resolution_clock::now();
	int32 kk = 0;
	
	while (true) {

		std::vector<Result> res;
		double minc = std::numeric_limits<double>::infinity();
		const auto start1 = std::chrono::high_resolution_clock::now();
		for (const auto& b : Boxes) {

			const FBox aabb = b->GetDefaultObject<APackerBox>()->get_aabb();
			//const auto start = std::chrono::high_resolution_clock::now();
			//Z_XY
			if constexpr (true) {
				const auto start5 = std::chrono::high_resolution_clock::now();
				auto ro = overlap(std::ceil(aabb.GetExtent().X), std::ceil(aabb.GetExtent().Y), std::ceil(aabb.GetExtent().Z));
				const std::chrono::duration<double> e5 = std::chrono::high_resolution_clock::now() - start5;
				s5 += e5.count();
				kk++;
				if (ro) {
					Result& r = ro.value();
					if (r.weight < minc) {
						minc = std::min(r.weight, minc);
						r.box = b;
						r.ext = FVector(std::ceil(aabb.GetExtent().X), std::ceil(aabb.GetExtent().Y), std::ceil(aabb.GetExtent().Z));
						r.perm = EAxisPerm::Z_XY_0;
						res.push_back(r);
					} //else continue;
				}
			}
			//Z_YX
			if constexpr (true) {
				auto ro = overlap(std::ceil(aabb.GetExtent().Y), std::ceil(aabb.GetExtent().X), std::ceil(aabb.GetExtent().Z));
				if (ro) {
					Result& r = ro.value();
					if (r.weight < minc) {
						minc = std::min(r.weight, minc);
						r.box = b;
						r.ext = FVector(std::ceil(aabb.GetExtent().Y), std::ceil(aabb.GetExtent().X), std::ceil(aabb.GetExtent().Z));
						r.perm = EAxisPerm::Z_XY_1;
						res.push_back(r);
					} //else continue;
				}
			}
			//Y_XZ
			if constexpr (true) {
				auto ro = overlap(std::ceil(aabb.GetExtent().X), std::ceil(aabb.GetExtent().Z), std::ceil(aabb.GetExtent().Y));
				if (ro) {
					Result& r = ro.value();
					if (r.weight < minc) {
						minc = std::min(r.weight, minc);
						r.box = b;
						r.ext = FVector(std::ceil(aabb.GetExtent().X), std::ceil(aabb.GetExtent().Z), std::ceil(aabb.GetExtent().Y));
						r.perm = EAxisPerm::Y_XZ_0;
						res.push_back(r);
					} //else continue;
				}
			}
			//Y_ZX
			if constexpr (true) {
				auto ro = overlap(std::ceil(aabb.GetExtent().Z), std::ceil(aabb.GetExtent().X), std::ceil(aabb.GetExtent().Y));
				if (ro) {
					Result& r = ro.value();
					if (r.weight < minc) {
						minc = std::min(r.weight, minc);
						r.box = b;
						r.ext = FVector(std::ceil(aabb.GetExtent().Z), std::ceil(aabb.GetExtent().X), std::ceil(aabb.GetExtent().Y));
						r.perm = EAxisPerm::Y_XZ_1;
						res.push_back(r);
					} //else continue;
				}
			}
			//X_YZ
			if constexpr (true) {
				auto ro = overlap(std::ceil(aabb.GetExtent().Y), std::ceil(aabb.GetExtent().Z), std::ceil(aabb.GetExtent().X));
				if (ro) {
					Result& r = ro.value();
					if (r.weight < minc) {
						minc = std::min(r.weight, minc);
						r.box = b;
						r.ext = FVector(std::ceil(aabb.GetExtent().Y), std::ceil(aabb.GetExtent().Z), std::ceil(aabb.GetExtent().X));
						r.perm = EAxisPerm::X_YZ_0;
						res.push_back(r);
					} //else continue;
				}
			}
			//X_ZY
			if constexpr(true){
				auto ro = overlap(std::ceil(aabb.GetExtent().Z), std::ceil(aabb.GetExtent().Y), std::ceil(aabb.GetExtent().X));
				if (ro) {
					Result& r = ro.value();
					if (r.weight < minc) {
						minc = std::min(r.weight, minc);
						r.box = b;
						r.ext = FVector(std::ceil(aabb.GetExtent().Z), std::ceil(aabb.GetExtent().Y), std::ceil(aabb.GetExtent().X));
						r.perm = EAxisPerm::X_YZ_1;
						res.push_back(r);
					} //else continue;
				}
			}

			//const std::chrono::duration<double> ee = std::chrono::high_resolution_clock::now() - start;
			//std::cout << "r: " << ee.count() << std::endl;

		}

		const std::chrono::duration<double> e1 = std::chrono::high_resolution_clock::now() - start1;
		s1 += e1.count();

		//if (k == 1) return;

		if (res.empty()) break;

		const auto start2 = std::chrono::high_resolution_clock::now();
		std::sort(res.begin(), res.end(), [](const Result& _e1, const Result& _e2) {
			return _e1.weight < _e2.weight;
		});
		const std::chrono::duration<double> e2 = std::chrono::high_resolution_clock::now() - start2;
		s2 += e2.count();

		const auto& r = res[0];
		APackerBox* box = r.box->GetDefaultObject<APackerBox>();

		const FBox tar = FBox(
			FVector(r.n0 - r.ext.X, r.n1 - r.ext.Y, r.h), 
			FVector(r.n0 + r.ext.X, r.n1 + r.ext.Y, r.h + 2*r.ext.Z));

		DrawDebugBox(world, tar.GetCenter(), tar.GetExtent(), FColor::Red, true);

		vol += tar.GetVolume();

		const auto start3 = std::chrono::high_resolution_clock::now();
		toSpawn.emplace_back(box->GetClass(), ::Detail::make_transform(
			r.perm,
			tar,
			box->get_aabb().GetCenter(),
			box->get_relative_location()
		));
		const std::chrono::duration<double> e3 = std::chrono::high_resolution_clock::now() - start3;
		s3 += e3.count();
		
		k++;

		//if (k > 50) break;

		//std::cout << "----" << std::endl;
		//std::cout << int32(tar.Min.X) / Tree::BUCKET_SIZE << ", " << int32(tar.Max.X) / Tree::BUCKET_SIZE  << std::endl;
		//std::cout << int32(tar.Min.Y) / Tree::BUCKET_SIZE << ", " << int32(tar.Max.Y) / Tree::BUCKET_SIZE << std::endl;
		
		std::fill(mm.begin(), mm.end(), false);
		for (int32 n0 = int32(tar.Min.X) / Tree::BUCKET_SIZE; n0 <= int32(tar.Max.X) / Tree::BUCKET_SIZE; ++n0) {
			for (int32 n1 = int32(tar.Min.Y) / Tree::BUCKET_SIZE; n1 <= int32(tar.Max.Y) / Tree::BUCKET_SIZE; ++n1) {
				const int32_t iid = n0 + n1 * bc;
				mm[iid] = true;
			}
		}

		for (int32 n0 = int32(tar.Min.X); n0 <= int32(tar.Max.X); ++n0) {
			for (int32 n1 = int32(tar.Min.Y); n1 <= int32(tar.Max.Y); ++n1) {
				const int32 i = n1 + n0 * Bounds;
				//if (i < 0 || i >= map.size()) {
				//	std::cout << "insert oob" << std::endl;
				//}
				map[i] = tar.Max.Z;
			}
		}

		//const auto start = std::chrono::high_resolution_clock::now();
		const auto start4 = std::chrono::high_resolution_clock::now();
		tree.recompute(mm);
		const std::chrono::duration<double> e4 = std::chrono::high_resolution_clock::now() - start4;
		s4 += e4.count();
		//const std::chrono::duration<double> ee = std::chrono::high_resolution_clock::now() - start;
		//std::cout << "r: " << ee.count() << std::endl;

		//std::stringstream ss;
		//ss << "map_" << k++ << ".png";
		//image_real<double, false>(Bounds, Bounds, map.data(), ss.str());
	}

	const std::chrono::duration<double> eet = std::chrono::high_resolution_clock::now() - startt;
	std::cout << "Spawned " << k << " boxes in " << eet.count() << "s" << std::endl;
	std::cout << "overlap: " << s1 << "s" << std::endl;
	std::cout << "sort: " << s2 << "s" << std::endl;
	std::cout << "transform: " << s3 << "s" << std::endl;
	std::cout << "rebuild: " << s4 << "s" << std::endl;
	std::cout << "#overlaps: " << kkk << std::endl;
	std::cout << "t/overlap: " << s1 / double(kkk) << "s" << std::endl;
	std::cout << "Volume: " << (1. / double(Bounds * Bounds * Height) * 100. * vol) << "%" << std::endl;

	FActorSpawnParameters params;
	params.bHideFromSceneOutliner = true;
	for (const auto& [c, t] : toSpawn) {
		world->SpawnActor<AActor>(c, t, params);
	}
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