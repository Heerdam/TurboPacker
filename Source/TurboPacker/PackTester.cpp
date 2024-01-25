
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

double APackerBox::get_weight() const {
	check(mesh);
	const auto aabb = get_aabb();
	return aabb.GetVolume();
}//APackerBox::get_weight

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
			APackerBox* box = world->SpawnActor<APackerBox>(c, t, params);
		}
		toSpawn.clear();
	}

}//APackTester::Tick

void APackTester::Clear() {
	UWorld* world = GetWorld();
	if (!world) return;

	{
		std::lock_guard<std::mutex> lock(*m);
		isPacking = false;
	}
	if (future != nullptr) {
		future->Wait();
	}
	toSpawn.clear();

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

				if (h != 0) continue;

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

		std::mutex mut;
		std::vector<::Detail::Result> res;
		std::atomic<double> minc = std::numeric_limits<double>::infinity();

		std::atomic<int32> c = 0;

		for (const auto& b : Boxes) {

			const FBox aabb = b->GetDefaultObject<APackerBox>()->get_aabb();
			c += 6;

			//Z_XY
			AsyncTask(ENamedThreads::AnyHiPriThreadHiPriTask, 
				[&, _b = b, _perm = EAxisPerm::Z_XY_0, 
				_ext0 = std::ceil(aabb.GetExtent().X), 
				_ext1 = std::ceil(aabb.GetExtent().Y), 
				_h = std::ceil(aabb.GetExtent().Z)]() -> void {	
					auto ro = overlap(_ext0, _ext1, _h);
					if (ro) {
						::Detail::Result& r = ro.value();
						if (r.weight < minc) {
							minc = std::min(r.weight, minc.load());
							r.box = _b;
							r.ext = FVector(_ext0, _ext1, _h);
							r.perm = _perm;
							std::lock_guard<std::mutex> l(mut);
							res.push_back(r);
						}
					}
					c--;
				}
			);

			//Z_YX
			AsyncTask(ENamedThreads::AnyHiPriThreadHiPriTask, 
				[&, _b = b, _perm = EAxisPerm::Z_XY_1,
					_ext0 = std::ceil(aabb.GetExtent().Y),
					_ext1 = std::ceil(aabb.GetExtent().X),
					_h = std::ceil(aabb.GetExtent().Z)]() -> void {
					auto ro = overlap(_ext0, _ext1, _h);
					if (ro) {
						::Detail::Result& r = ro.value();
						if (r.weight < minc) {
							minc = std::min(r.weight, minc.load());
							r.box = _b;
							r.ext = FVector(_ext0, _ext1, _h);
							r.perm = _perm;
							std::lock_guard<std::mutex> l(mut);
							res.push_back(r);
						}
					}
					c--;
				}
				);

			//Y_XZ
			AsyncTask(ENamedThreads::AnyHiPriThreadHiPriTask, 
				[&, _b = b, _perm = EAxisPerm::Y_XZ_0,
					_ext0 = std::ceil(aabb.GetExtent().X),
					_ext1 = std::ceil(aabb.GetExtent().Z),
					_h = std::ceil(aabb.GetExtent().Y)]() -> void {
					auto ro = overlap(_ext0, _ext1, _h);
					if (ro) {
						::Detail::Result& r = ro.value();
						if (r.weight < minc) {
							minc = std::min(r.weight, minc.load());
							r.box = _b;
							r.ext = FVector(_ext0, _ext1, _h);
							r.perm = _perm;
							std::lock_guard<std::mutex> l(mut);
							res.push_back(r);
						}
					}
					c--;
				}
				);

			//Y_ZX
			AsyncTask(ENamedThreads::AnyHiPriThreadHiPriTask, 
				[&, _b = b, _perm = EAxisPerm::Y_XZ_1,
					_ext0 = std::ceil(aabb.GetExtent().Z),
					_ext1 = std::ceil(aabb.GetExtent().X),
					_h = std::ceil(aabb.GetExtent().Y)]() -> void {
					auto ro = overlap(_ext0, _ext1, _h);
					if (ro) {
						::Detail::Result& r = ro.value();
						if (r.weight < minc) {
							minc = std::min(r.weight, minc.load());
							r.box = _b;
							r.ext = FVector(_ext0, _ext1, _h);
							r.perm = _perm;
							std::lock_guard<std::mutex> l(mut);
							res.push_back(r);
						}
					}
					c--;
				}
				);

			//X_YZ
			AsyncTask(ENamedThreads::AnyHiPriThreadHiPriTask, 
				[&, _b = b, _perm = EAxisPerm::X_YZ_0,
					_ext0 = std::ceil(aabb.GetExtent().Y),
					_ext1 = std::ceil(aabb.GetExtent().Z),
					_h = std::ceil(aabb.GetExtent().X)]() -> void {
					auto ro = overlap(_ext0, _ext1, _h);
					if (ro) {
						::Detail::Result& r = ro.value();
						if (r.weight < minc) {
							minc = std::min(r.weight, minc.load());
							r.box = _b;
							r.ext = FVector(_ext0, _ext1, _h);
							r.perm = _perm;
							std::lock_guard<std::mutex> l(mut);
							res.push_back(r);
						}
					}
					c--;
				}
				);

			//X_ZY
			AsyncTask(ENamedThreads::AnyHiPriThreadHiPriTask, 
				[&, _b = b, _perm = EAxisPerm::X_YZ_1,
					_ext0 = std::ceil(aabb.GetExtent().Z),
					_ext1 = std::ceil(aabb.GetExtent().Y),
					_h = std::ceil(aabb.GetExtent().X)]() -> void {
					auto ro = overlap(_ext0, _ext1, _h);
					if (ro) {
						::Detail::Result& r = ro.value();
						if (r.weight < minc) {
							minc = std::min(r.weight, minc.load());
							r.box = _b;
							r.ext = FVector(_ext0, _ext1, _h);
							r.perm = _perm;
							std::lock_guard<std::mutex> l(mut);
							res.push_back(r);
						}
					}
					c--;
				}
				);

		}

		while(c != 0){}

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
	std::cout << "start" << std::endl;
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

// Naive:
//  Spawned 162 boxes in 57.0335s - 0.352s /box
// MQT:
//  Spawned 162 boxes in 15.7559s - 0.09725s/box
//  Spawned 648 boxes in 160.478s - 0.247s/box

//---------------------------------------------------

AOnlinePacker::AOnlinePacker() {

	m = std::make_unique<std::mutex>();

	PrimaryActorTick.bCanEverTick = true;
	PrimaryActorTick.TickInterval = 1.f;

}//AOnlinePacker::AOnlinePacker

void AOnlinePacker::Tick(float _delta) {

	if (isPacking) {
		std::lock_guard<std::mutex> lock(*m);

		UWorld* world = GetWorld();
		FActorSpawnParameters params;
		params.bHideFromSceneOutliner = true;
		for (const auto& [c, t] : toSpawn) {
			APackerBox* box = world->SpawnActor<APackerBox>(c, t, params);
			vol += box->get_aabb().GetVolume();
		}
		toSpawn.clear();
	}

}//AOnlinePacker::Tick

void AOnlinePacker::Clear() {
	UWorld* world = GetWorld();
	if (!world) return;

	{
		std::lock_guard<std::mutex> lock(*m);
		isPacking = false;
	}
	if (future != nullptr) {
		future->Wait();
	}
	toSpawn.clear();
	q = std::queue<APackerBox*>();
	vol = 0.;

	FlushPersistentDebugLines(world);

	TArray<AActor*> tod;
	UGameplayStatics::GetAllActorsOfClass(world, APackerBox::StaticClass(), tod);
	for (AActor* a : tod)
		world->DestroyActor(a);
}//AOnlinePacker::Clear

void AOnlinePacker::pack_impl() {

	using namespace MQT2;

	UOnlinePackerConfig* conf = Config->GetDefaultObject<UOnlinePackerConfig>();

	const auto overlap = [&](const int32 _ext0, const int32 _ext1, const int32 _h)->std::optional<::Detail::Result> {

		using namespace ::Detail;

		for (int32 n0 = _ext0; n0 < conf->Bounds - _ext0; ++n0) {
			for (int32 n1 = _ext1; n1 < conf->Bounds - _ext1; ++n1) {

				const int32 i = n1 + n0 * conf->Bounds;
				if (map[i] + 2 * _h >= conf->Height) continue;

				const auto [l, m, h] = tree->check_overlap(
					Vec2{ int32_t(n0 - _ext0), int32_t(n1 - _ext1) },
					Vec2{ int32_t(n0 + _ext0), int32_t(n1 + _ext1) },
					map[i]);

				if (h != 0 ) continue;

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

	const int32_t bc = (conf->Bounds / Tree::BUCKET_SIZE);
	std::vector<bool> mm;
	mm.resize(bc * bc);
	std::fill(mm.begin(), mm.end(), true);

	while (!q.empty()) {

		std::mutex mut;
		std::vector<::Detail::Result> res;
		std::atomic<double> minc = std::numeric_limits<double>::infinity();

		std::atomic<int32> c = 0;

		APackerBox* next = q.front();
		q.pop();

		const auto b = next->GetClass();

		const FBox aabb = next->get_aabb();
		c += 6;

		//Z_XY
		AsyncTask(ENamedThreads::AnyHiPriThreadHiPriTask, 
			[&, _b = b, _perm = EAxisPerm::Z_XY_0,
				_ext0 = std::ceil(aabb.GetExtent().X),
				_ext1 = std::ceil(aabb.GetExtent().Y),
				_h = std::ceil(aabb.GetExtent().Z)]() -> void {
				auto ro = overlap(_ext0, _ext1, _h);
				if (ro) {
					::Detail::Result& r = ro.value();
					if (r.weight < minc) {
						minc = std::min(r.weight, minc.load());
						r.box = _b;
						r.ext = FVector(_ext0, _ext1, _h);
						r.perm = _perm;
						std::lock_guard<std::mutex> l(mut);
						res.push_back(r);
					}
				}
				c--;
			}
			);

		//Z_YX
		AsyncTask(ENamedThreads::AnyHiPriThreadHiPriTask, 
			[&, _b = b, _perm = EAxisPerm::Z_XY_1,
				_ext0 = std::ceil(aabb.GetExtent().Y),
				_ext1 = std::ceil(aabb.GetExtent().X),
				_h = std::ceil(aabb.GetExtent().Z)]() -> void {
				auto ro = overlap(_ext0, _ext1, _h);
				if (ro) {
					::Detail::Result& r = ro.value();
					if (r.weight < minc) {
						minc = std::min(r.weight, minc.load());
						r.box = _b;
						r.ext = FVector(_ext0, _ext1, _h);
						r.perm = _perm;
						std::lock_guard<std::mutex> l(mut);
						res.push_back(r);
					}
				}
				c--;
			}
			);

		//Y_XZ
		AsyncTask(ENamedThreads::AnyHiPriThreadHiPriTask, 
			[&, _b = b, _perm = EAxisPerm::Y_XZ_0,
				_ext0 = std::ceil(aabb.GetExtent().X),
				_ext1 = std::ceil(aabb.GetExtent().Z),
				_h = std::ceil(aabb.GetExtent().Y)]() -> void {
				auto ro = overlap(_ext0, _ext1, _h);
				if (ro) {
					::Detail::Result& r = ro.value();
					if (r.weight < minc) {
						minc = std::min(r.weight, minc.load());
						r.box = _b;
						r.ext = FVector(_ext0, _ext1, _h);
						r.perm = _perm;
						std::lock_guard<std::mutex> l(mut);
						res.push_back(r);
					}
				}
				c--;
			}
			);

		//Y_ZX
		AsyncTask(ENamedThreads::AnyHiPriThreadHiPriTask, 
			[&, _b = b, _perm = EAxisPerm::Y_XZ_1,
				_ext0 = std::ceil(aabb.GetExtent().Z),
				_ext1 = std::ceil(aabb.GetExtent().X),
				_h = std::ceil(aabb.GetExtent().Y)]() -> void {
				auto ro = overlap(_ext0, _ext1, _h);
				if (ro) {
					::Detail::Result& r = ro.value();
					if (r.weight < minc) {
						minc = std::min(r.weight, minc.load());
						r.box = _b;
						r.ext = FVector(_ext0, _ext1, _h);
						r.perm = _perm;
						std::lock_guard<std::mutex> l(mut);
						res.push_back(r);
					}
				}
				c--;
			}
			);

		//X_YZ
		AsyncTask(ENamedThreads::AnyHiPriThreadHiPriTask, 
			[&, _b = b, _perm = EAxisPerm::X_YZ_0,
				_ext0 = std::ceil(aabb.GetExtent().Y),
				_ext1 = std::ceil(aabb.GetExtent().Z),
				_h = std::ceil(aabb.GetExtent().X)]() -> void {
				auto ro = overlap(_ext0, _ext1, _h);
				if (ro) {
					::Detail::Result& r = ro.value();
					if (r.weight < minc) {
						minc = std::min(r.weight, minc.load());
						r.box = _b;
						r.ext = FVector(_ext0, _ext1, _h);
						r.perm = _perm;
						std::lock_guard<std::mutex> l(mut);
						res.push_back(r);
					}
				}
				c--;
			}
			);

		//X_ZY
		AsyncTask(ENamedThreads::AnyHiPriThreadHiPriTask, 
			[&, _b = b, _perm = EAxisPerm::X_YZ_1,
				_ext0 = std::ceil(aabb.GetExtent().Z),
				_ext1 = std::ceil(aabb.GetExtent().Y),
				_h = std::ceil(aabb.GetExtent().X)]() -> void {
				auto ro = overlap(_ext0, _ext1, _h);
				if (ro) {
					::Detail::Result& r = ro.value();
					if (r.weight < minc) {
						minc = std::min(r.weight, minc.load());
						r.box = _b;
						r.ext = FVector(_ext0, _ext1, _h);
						r.perm = _perm;
						std::lock_guard<std::mutex> l(mut);
						res.push_back(r);
					}
				}
				c--;
			}
			);

		while (c != 0) {}

		if (!isPacking) break;
		if (res.empty()) continue;

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
				const int32 i = n1 + n0 * conf->Bounds;
				map[i] = tar.Max.Z;
			}
		}

		tree->recompute(mm);

	}

}//AOnlinePacker::pack_impl

void AOnlinePacker::Pack() {
	UWorld* world = GetWorld();
	if (!world) return;

	if (Config == nullptr) {
		std::cout << "Error: no Config set!" << std::endl;
		return;
	}

	Clear();

	using namespace MQT2;

	UOnlinePackerConfig* conf = Config->GetDefaultObject<UOnlinePackerConfig>();

	map.resize(conf->Bounds * conf->Bounds);
	std::fill(map.begin(), map.end(), 0.);

	tree = std::make_unique<Tree>(map, conf->Bounds);

	DrawDebugBox(world,
		FVector(conf->Bounds * 0.5, conf->Bounds * 0.5, conf->Height * 0.5),
		FVector(conf->Bounds * 0.5, conf->Bounds * 0.5, conf->Height * 0.5), FColor::Blue, true);

	//---------------

	std::vector<APackerBox*> qq;
	for (const auto& p : conf->Boxes) {
		for (int32 i = 0; i < p.Count; ++i)
			qq.push_back(p.Type->GetDefaultObject<APackerBox>());
	}

	std::sort(qq.begin(), qq.end(), [](APackerBox* _o1, APackerBox* _o2) {
		return _o1->get_weight() > _o2->get_weight();

	});

	for (const auto& ptr : qq)
		q.push(ptr);

	//---------------

	isPacking = true;
	std::cout << "start" << std::endl;
	start = std::chrono::high_resolution_clock::now();
	future = std::make_unique<TFuture<bool>>(AsyncThread([this]() -> bool {
		pack_impl();
		const std::chrono::duration<double> ee = std::chrono::high_resolution_clock::now() - start;
		std::cout << "done (" << ee.count() << "s)" << std::endl;

		UOnlinePackerConfig* conf = Config->GetDefaultObject<UOnlinePackerConfig>();
		std::cout << 1. / double(conf->Bounds * conf->Bounds * conf->Height) * vol << std::endl;
		return true;
		})
	);

}//AOnlinePacker::Pack