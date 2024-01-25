
#include "PackTester.h"


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

APackTester::APackTester() {

	m = std::make_unique<std::mutex>();

	PrimaryActorTick.bCanEverTick = true;
	PrimaryActorTick.TickInterval = 1.f;

}//APackTester::APackTester

void APackTester::Tick(float _delta) {

	if(isPacking){
		std::lock_guard<std::mutex> lock(*m);
		
		UWorld* world = GetWorld();
		FActorSpawnParameters params;
		params.bHideFromSceneOutliner = true;
		for (const auto& [c, t] : toSpawn) {
			world->SpawnActor<AActor>(c, t, params);
		}
		toSpawn.clear();
	}

}//APackTester::Tick

void APackTester::Clear() {
	UWorld* world = GetWorld();
	if (!world) return;

	{
		std::lock_guard<std::mutex> lock(*m);
		if (future != nullptr) {
			future->Wait();
		}
		isPacking = false;
		toSpawn.clear();

	}

	FlushPersistentDebugLines(world);
	
	TArray<AActor*> tod;
	UGameplayStatics::GetAllActorsOfClass(world, APackerBox::StaticClass(), tod);
	for (AActor* a : tod)
		world->DestroyActor(a);
}

void APackTester::pack_impl() {

	using namespace MQT2;

	const auto overlap = [&](const int32 _ext0, const int32 _ext1, const int32 _h)->std::optional<::Detail::Result> {

		using namespace ::Detail;

		for (int32 n0 = _ext0; n0 < Bounds - _ext0; ++n0) {
			for (int32 n1 = _ext1; n1 < Bounds - _ext1; ++n1) {

				const int32 i = n1 + n0 * Bounds;
				if (map[i] + 2 * _h >= Height) continue;

				const auto [l, m, h] = tree->check_overlap(
					Vec2{ int32_t(n0 - _ext0), int32_t(n1 - _ext1) },
					Vec2{ int32_t(n0 + _ext0), int32_t(n1 + _ext1) },
					map[i]);

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
	mm.resize(bc * bc);
	std::fill(mm.begin(), mm.end(), true);

	while (true) {

		std::vector<::Detail::Result> res;
		double minc = std::numeric_limits<double>::infinity();

		for (const auto& b : Boxes) {

			const FBox aabb = b->GetDefaultObject<APackerBox>()->get_aabb();

			//Z_XY
			if constexpr (true) {
				auto ro = overlap(std::ceil(aabb.GetExtent().X), std::ceil(aabb.GetExtent().Y), std::ceil(aabb.GetExtent().Z));
				if (ro) {
					::Detail::Result& r = ro.value();
					if (r.weight < minc) {
						minc = std::min(r.weight, minc);
						r.box = b;
						r.ext = FVector(std::ceil(aabb.GetExtent().X), std::ceil(aabb.GetExtent().Y), std::ceil(aabb.GetExtent().Z));
						r.perm = EAxisPerm::Z_XY_0;
						res.push_back(r);
					}
				}
			}
			//Z_YX
			if constexpr (true) {
				auto ro = overlap(std::ceil(aabb.GetExtent().Y), std::ceil(aabb.GetExtent().X), std::ceil(aabb.GetExtent().Z));
				if (ro) {
					::Detail::Result& r = ro.value();
					if (r.weight < minc) {
						minc = std::min(r.weight, minc);
						r.box = b;
						r.ext = FVector(std::ceil(aabb.GetExtent().Y), std::ceil(aabb.GetExtent().X), std::ceil(aabb.GetExtent().Z));
						r.perm = EAxisPerm::Z_XY_1;
						res.push_back(r);
					}
				}
			}
			//Y_XZ
			if constexpr (true) {
				auto ro = overlap(std::ceil(aabb.GetExtent().X), std::ceil(aabb.GetExtent().Z), std::ceil(aabb.GetExtent().Y));
				if (ro) {
					::Detail::Result& r = ro.value();
					if (r.weight < minc) {
						minc = std::min(r.weight, minc);
						r.box = b;
						r.ext = FVector(std::ceil(aabb.GetExtent().X), std::ceil(aabb.GetExtent().Z), std::ceil(aabb.GetExtent().Y));
						r.perm = EAxisPerm::Y_XZ_0;
						res.push_back(r);
					}
				}
			}
			//Y_ZX
			if constexpr (true) {
				auto ro = overlap(std::ceil(aabb.GetExtent().Z), std::ceil(aabb.GetExtent().X), std::ceil(aabb.GetExtent().Y));
				if (ro) {
					::Detail::Result& r = ro.value();
					if (r.weight < minc) {
						minc = std::min(r.weight, minc);
						r.box = b;
						r.ext = FVector(std::ceil(aabb.GetExtent().Z), std::ceil(aabb.GetExtent().X), std::ceil(aabb.GetExtent().Y));
						r.perm = EAxisPerm::Y_XZ_1;
						res.push_back(r);
					}
				}
			}
			//X_YZ
			if constexpr (true) {
				auto ro = overlap(std::ceil(aabb.GetExtent().Y), std::ceil(aabb.GetExtent().Z), std::ceil(aabb.GetExtent().X));
				if (ro) {
					::Detail::Result& r = ro.value();
					if (r.weight < minc) {
						minc = std::min(r.weight, minc);
						r.box = b;
						r.ext = FVector(std::ceil(aabb.GetExtent().Y), std::ceil(aabb.GetExtent().Z), std::ceil(aabb.GetExtent().X));
						r.perm = EAxisPerm::X_YZ_0;
						res.push_back(r);
					}
				}
			}
			//X_ZY
			if constexpr (true) {
				auto ro = overlap(std::ceil(aabb.GetExtent().Z), std::ceil(aabb.GetExtent().Y), std::ceil(aabb.GetExtent().X));
				if (ro) {
					::Detail::Result& r = ro.value();
					if (r.weight < minc) {
						minc = std::min(r.weight, minc);
						r.box = b;
						r.ext = FVector(std::ceil(aabb.GetExtent().Z), std::ceil(aabb.GetExtent().Y), std::ceil(aabb.GetExtent().X));
						r.perm = EAxisPerm::X_YZ_1;
						res.push_back(r);
					}
				}
			}

		}


		if (res.empty() || !isPacking) break;

		std::sort(res.begin(), res.end(), [](const ::Detail::Result& _e1, const ::Detail::Result& _e2) {
			return _e1.weight < _e2.weight;
		});

		const auto& r = res[0];
		APackerBox* box = r.box->GetDefaultObject<APackerBox>();

		const FBox tar = FBox(
			FVector(r.n0 - r.ext.X, r.n1 - r.ext.Y, r.h),
			FVector(r.n0 + r.ext.X, r.n1 + r.ext.Y, r.h + 2 * r.ext.Z));

		{
			std::lock_guard<std::mutex> lock(*m);
			toSpawn.emplace_back(box->GetClass(), ::Detail::make_transform(
				r.perm,
				tar,
				box->get_aabb().GetCenter(),
				box->get_relative_location()
			));
		}

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
				map[i] = tar.Max.Z;
			}
		}

		tree->recompute(mm);

	}

}//APackTester::pack_impl

void APackTester::Pack() {
	UWorld* world = GetWorld();
	if (!world) return;

	Clear();

	using namespace MQT2;

	map.resize(Bounds * Bounds);
	std::fill(map.begin(), map.end(), 0.);

	tree = std::make_unique<Tree>(map, Bounds);

	DrawDebugBox(world, 
		FVector(Bounds * 0.5, Bounds * 0.5, Height * 0.5),
		FVector(Bounds * 0.5, Bounds * 0.5, Height * 0.5), FColor::Blue, true);

	//---------------

	isPacking = true;

	future = std::make_unique<TFuture<bool>>(AsyncThread([this]() -> bool {
		pack_impl();
		std::cout << "done" << std::endl;
		return true;
		})
	);

	/*const std::chrono::duration<double> eet = std::chrono::high_resolution_clock::now() - startt;
	std::cout << "Spawned " << k << " boxes in " << eet.count() << "s" << std::endl;
	std::cout << "overlap: " << s1 << "s" << std::endl;
	std::cout << "sort: " << s2 << "s" << std::endl;
	std::cout << "transform: " << s3 << "s" << std::endl;
	std::cout << "rebuild: " << s4 << "s" << std::endl;
	std::cout << "#overlaps: " << kkk << std::endl;
	std::cout << "t/overlap: " << s1 / double(kkk) << "s" << std::endl;
	std::cout << "Volume: " << (1. / double(Bounds * Bounds * Height) * 100. * vol) << "%" << std::endl;*/


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