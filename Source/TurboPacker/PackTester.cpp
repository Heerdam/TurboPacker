
#include "PackTester.h"

#include "Materials/MaterialInstanceDynamic.h"
#include "Kismet/GameplayStatics.h"
#include "EnhancedActionKeyMapping.h"
#include "EnhancedInputComponent.h"
#include "EnhancedInputSubsystems.h"
#include "InputMappingContext.h"

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

FBox APackerBox::get_aabb(const FVector) const {
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
	const auto aabb = get_aabb(FVector());
	return aabb.GetVolume();
}//APackerBox::get_weight

//----------------------------------------

ARandomBox::ARandomBox() {}//ARandomBox::ARandomBox

void ARandomBox::set_to_size(const FVector& _new_size, const double _min_size, const double _max_size, bool _isCube) {
	check(mesh);
	mesh->SetWorldScale3D(_new_size / 100.);
	UMaterialInstanceDynamic* mi = UMaterialInstanceDynamic::Create(mesh->GetMaterial(0), this);
	const double temp = 1. - 1. / (_max_size - _min_size) * ((_new_size.X * _new_size.Y * _new_size.Z) - _min_size);
	//std::cout << _min_size << ", " << _max_size << ", " << (_new_size.X * _new_size.Y * _new_size.Z) << std::endl;
	//std::cout << temp << std::endl;
	//const FLinearColor c1 = FLinearColor::Blue;
	//const FLinearColor c2 = FLinearColor::Green;
	//const FLinearColor c3 = FLinearColor::Red;
	//const FLinearColor res = temp <= 0.5 ? FLinearColor::LerpUsingHSV(c1, c2, temp * 2.) : FLinearColor::LerpUsingHSV(c2, c3, (temp - 0.5) * 2.);
	mi->SetVectorParameterValue(TEXT("Color"), FColor::MakeRedToGreenColorFromScalar(temp));
	mesh->SetMaterial(0, mi);

}//ARandomBox::set_to_size

FBox ARandomBox::get_aabb(const FVector _ext) const {
	return { FVector(-_ext.X, -_ext.Y, -_ext.Z), FVector(_ext.X, _ext.Y, _ext.Z) };
}

//----------------------------------------
//----------------------------------------
//----------------------------------------

void AObserverController::BeginPlay() {
	Super::BeginPlay();

	UWorld* world = GetWorld();
	if (world) {
		TArray<AActor*> p;
		UGameplayStatics::GetAllActorsOfClass(world, APacker::StaticClass(), p);
		if (!p.IsEmpty())
			Packer = Cast<APacker>(p[0]);
	}
}//AObserverController::BeginPlay

void AObserverController::OnPossess(APawn* aPawn) {
	Super::OnPossess(aPawn);

	if (UEnhancedInputLocalPlayerSubsystem* Subsystem = ULocalPlayer::GetSubsystem<UEnhancedInputLocalPlayerSubsystem>(GetLocalPlayer())) {
		Subsystem->AddMappingContext(IA_context, 0);
	}

}//AObserverController::OnPossess

void AObserverController::SetupInputComponent() {
	Super::SetupInputComponent();

	UEnhancedInputComponent* input = CastChecked<UEnhancedInputComponent>(InputComponent);
	if (input) {
		//Pack
		input->BindAction(IA_context->GetMapping(0).Action, ETriggerEvent::Triggered, this, &AObserverController::Pack);
		
		//Clear
		input->BindAction(IA_context->GetMapping(1).Action, ETriggerEvent::Triggered, this, &AObserverController::Clear);
	}
}//AObserverController::SetupInputComponent

void AObserverController::Pack(const FInputActionValue& Value) {
	std::cout << "ding" << std::endl;
	if (Packer) {
		Packer->Pack();
	}
}// AObserverController::Pack

void AObserverController::Clear(const FInputActionValue& Value) {
	if (Packer) {
		Packer->Clear();
	}
}//AObserverController::Clear

//---------------------------------------------------
//---------------------------------------------------
//---------------------------------------------------

APacker::APacker() {

	m_ = std::make_unique<std::mutex>();

	PrimaryActorTick.bCanEverTick = true;
	PrimaryActorTick.TickInterval = 1.f;

}//APacker::APacker

double APacker::get_pack_percent() {
	UPackerConfig* conf = Config->GetDefaultObject<UPackerConfig>();
	if (conf) return 1. / double(conf->Bounds.X * conf->Bounds.Y * conf->Height) * vol_;
	else return 0.;
}//APacker::get_pack_percent

double APacker::get_time() {
	if (isPacking_) {
		const std::chrono::duration<double> ee = std::chrono::high_resolution_clock::now() - start_;
		return ee.count();
	} else return last_time_;
}//APacker::get_pack_percent

void APacker::Tick(float _delta) {
	UPackerConfig* conf = Config->GetDefaultObject<UPackerConfig>();
	UWorld* world = GetWorld();

	std::lock_guard<std::mutex> lock(*m_);
	for (const auto& res : toSpawn_) {

		FActorSpawnParameters params;
		params.bHideFromSceneOutliner = true;
		params.CustomPreSpawnInitalization = [&](AActor* _a) {
		ARandomBox* b = Cast<ARandomBox>(_a);
			b->set_to_size(2 * res.ext_org, conf->MinBoxVolume, conf->MaxBoxVolume, conf->CubeRandomBoxes);
		};

		if (res.isRandomBox) {
			ARandomBox* b = world->SpawnActor<ARandomBox>(res.box, res.trans, params);
				
		} else {
			world->SpawnActor<APackerBox>(res.box, res.trans, params);	
		}
		const FVector size = 2 * res.ext;
		vol_ += size.X * size.Y * size.Z;
		bcc_++;
	}
	toSpawn_.clear();

}//APacker::Tick

void APacker::Clear() {
	UWorld* world = GetWorld();
	if (!world) return;
	Stop();
	toSpawn_.clear();
	vol_ = 0.;
	bcc_ = 0;
	mcc_ = 0;

	FlushPersistentDebugLines(world);

	TArray<AActor*> tod;
	UGameplayStatics::GetAllActorsOfClass(world, APackerBox::StaticClass(), tod);
	for (AActor* a : tod)
		world->DestroyActor(a);
}//APacker::Clear

void APacker::Stop() {
	std::lock_guard<std::mutex> lock(*m_);
	isPacking_ = false;
	if (future_ != nullptr)
		future_->Wait();
}//APacker::Stop

std::vector<::Detail::Result> APacker::overlap_impl(
	const int32 _ext0, 
	const int32 _ext1, 
	const int32 _h
) {

	using namespace MQT2;
	using namespace ::Detail;

	UPackerConfig* conf = Config->GetDefaultObject<UPackerConfig>();
	const int32 ee0 = conf->Bounds.Y + 2;
	const int32 ee1 = conf->Bounds.X + 2;

	std::vector<::Detail::Result> res;
	for (int32 n0 = _ext0 + 1; n0 < ee0 - _ext0 - 1; ++n0) {
		for (int32 n1 = _ext1 + 1; n1 < ee1 - _ext1 - 1; ++n1) {

			const int32 i = n1 + n0 * N_;
			if (int32(map_[i]) + 2 * _h >= conf->Height) continue;

			const auto [l1, m1, h1] = tree_->check_overlap(
				Vec2{ int32_t(n0 - _ext0), int32_t(n1 - _ext1) },
				Vec2{ int32_t(n0 + _ext0), int32_t(n1 + _ext1) },
				map_[i]);

			const auto [l2, m2, h2] = tree_->check_border_overlap(
				Vec2{ int32_t(n0 - _ext0) - 1, int32_t(n1 - _ext1) + 1 },
				Vec2{ int32_t(n0 + _ext0) - 1, int32_t(n1 + _ext1) + 1 },
				map_[i]);

			if (conf->AllowOverlap) {
				if (h1 != 0) continue;
			} else {
				if (h1 != 0 || l1 != 0) continue;
			}

			Result out;
			out.n0 = n0;
			out.n1 = n1;
			out.l = l1;
			out.b_l = l2;
			out.b_m = m2;
			out.b_h = h2;
			out.h = map_[i];

			res.push_back(out);
		}
	}
	return res;

}//APacker::overlap_impl

void APacker::dispatch_impl(
	std::vector<::Detail::Result>& _res, 
	std::mutex& _m,
	std::atomic<double>& _minc,
	std::atomic<int32>& _c,
	const int32 _n0, const int32 _n1, 
	const FVector& _ext_org,
	const int32 _h, EAxisPerm _perm,
	const TSubclassOf<APackerBox>& _b,
	const std::function<double(const Detail::Result&)>& _cost) {

	UPackerConfig* conf = Config->GetDefaultObject<UPackerConfig>();

	auto tr = [&,
		_b = _b,
		_perm = _perm,
		_ext0 = _n0,
		_ext1 = _n1,
		_ext_org = _ext_org,
		_h = _h, cost = _cost] () -> void {
		UPackerConfig* conf = Config->GetDefaultObject<UPackerConfig>();
		auto ro = overlap_impl(_ext0, _ext1, _h);
		if (!ro.empty()) {
			for (auto& r : ro) {
				r.box = _b;
				r.ext = FVector(_ext0, _ext1, _h);
				r.ext_org = _ext_org;
				r.perm = _perm;
				r.isRandomBox = conf->BoxType == BoxGenerationType::RANDOM;
				r.weight = cost(r);
				if (r.weight < _minc) {
					_minc = std::min(r.weight, _minc.load());
					std::lock_guard<std::mutex> l(_m);
					_res.push_back(r);
				}
			}
		}
		_c--;
	};

	if (conf->MultiThreading) AsyncTask(ENamedThreads::AnyHiPriThreadHiPriTask, std::move(tr));
	else tr();

}//APacker::dispatch_impl

void APacker::pack_impl() {

	using namespace MQT2;

	UPackerConfig* conf = Config->GetDefaultObject<UPackerConfig>();

	const int32 ee0 = conf->Bounds.Y + 2;
	const int32 ee1 = conf->Bounds.X + 2;

	N_ = Util::get_power_of_2(std::max(ee0, ee1), Tree::BUCKET_SIZE);
	const bool useRandomBox = conf->BoxType == BoxGenerationType::RANDOM;

	map_.resize(N_ * N_);
	std::fill(map_.begin(), map_.end(), 0.);

	for (int32 n0 = 0; n0 < ee0; ++n0) {
		const int32 i1 = n0 * N_;
		const int32 i2 = (ee1 - 1) + n0 * N_;
		map_[i1] = conf->Height;
		map_[i2] = conf->Height;
	}

	for (int32 n1 = 0; n1 < ee1; ++n1) {
		const int32 i1 = n1;
		const int32 i2 = n1 + (ee0 - 1) * N_;
		map_[i1] = conf->Height;
		map_[i2] = conf->Height;
	}

	tree_ = std::make_unique<Tree>(map_, N_);

	DrawDebugBox(GetWorld(),
		FVector(ee0 * 0.5, ee0 * 0.5, conf->Height * 0.5),
		FVector(ee1 * 0.5, ee1 * 0.5, conf->Height * 0.5), FColor::Blue, true);

	//--------------------------

	const int32_t bc = (N_ / Tree::BUCKET_SIZE);
	std::vector<bool> mm;
	mm.resize(bc * bc);
	std::fill(mm.begin(), mm.end(), true);

	int32 et = 0;

	std::random_device rd;
	conf->Seed = conf->UseRandomSeed ? rd() : conf->Seed;
	Rand g(conf->Seed);
	DistD dist = DistD(conf->MinBoxVolume, conf->MaxBoxVolume);

	const auto co = [&](const ::Detail::Result& _r) {
		UPackerConfig* conf = Config->GetDefaultObject<UPackerConfig>();

		switch (conf->CostFunction) {
			case CostFunction::SIMPLE:
				return std::pow(double(_r.h), 3) + std::pow(_r.n0 + _r.n1, 2);
			case CostFunction::SIMPLE_HEIGHT:
				return std::pow(_r.n0 + _r.n1, 2) - std::pow(_r.b_h, 2);
		}
		return 0.;

		//return std::pow(double(_r.h), 3) + std::pow(_r.n0 + _r.n1, 2);
		//return std::pow(_r.n0 + _r.n1, 2) - std::pow(_r.b_h, 2);
		//return std::pow(double(_r.h), 3) + std::pow(std::abs((_r.n0 - conf->Bounds / 2)), 2) + std::pow(std::abs((_r.n1 - conf->Bounds / 2)), 2);
		//return std::pow(std::abs((_r.n0 - conf->Bounds / 2)), 2) + std::pow(std::abs((_r.n1 - conf->Bounds / 2)), 2);

		//const double mw = double(_r.ext.X * _r.ext.Y * _r.ext.Z);
		//const double temp = 1. / (conf->MaxBoxSize - conf->MinBoxSize) * (mw - conf->MinBoxSize);

		//if (temp < 0.5) {
		//	return std::pow(_r.n0 + _r.n1, 2);
		//} else {
		//	return std::pow((conf->Bounds - _r.n0) + (conf->Bounds - _r.n1), 2);
		//}
	};

	const auto rbox = [&](const double _vol, const bool _c) -> FVector {
		const double s = std::pow(_vol, 1. / 3.);
		if (_c) return FVector(s) * 0.5;
		const double ss = 3 * s;
		const double d1 = DistD(0.2 * ss, ss)(g);
		const double d2 = DistD(0.2 * (ss - d1), ss - d1)(g);
		const double d3 = ss - d1 - d2;
		return { d1, d2, d3 };
	};

	const auto b = conf->Box;
	const auto box = conf->Box->GetDefaultObject<ARandomBox>();

	std::vector<FVector> next_set(conf->SetSize);
	std::queue<FVector> next_q;

	if (conf->BoxType == BoxGenerationType::LIST) {
		UBoxList* list = conf->List[conf->ListIndex]->GetDefaultObject<UBoxList>();
		std::vector<FVector> tmp;
		for (const FPPair& p : list->List) {
			for (int32 i = 0; i < p.Count; ++i)
				tmp.push_back(p.Size);
		}
		if (list->ShuffleBoxes)
			std::shuffle(tmp.begin(), tmp.end(), g);
		for (const auto& f : tmp)
			next_q.push(f);
	}

	while (true) {

		if (conf->BoxType == BoxGenerationType::LIST && next_q.empty()) break;

		std::mutex mut;
		std::vector<::Detail::Result> res;
		std::atomic<double> minc = std::numeric_limits<double>::infinity();

		std::atomic<int32> c = 0;

		next_set.clear();
		
		switch (conf->BoxType) {
			case BoxGenerationType::RANDOM:
			{
				for (int32 i = 0; i < conf->SetSize; ++i)
					next_set.push_back(rbox(dist(g), conf->CubeRandomBoxes));
			}
			break;
			case BoxGenerationType::LIST:
			{
				for (int32 i = 0; i < conf->SetSize; ++i) {
					if (next_q.empty()) break;
					next_set.push_back(next_q.front());
					next_q.pop();
				}
			}
			break;
		}
		
		for (const FVector& nextSize : next_set) {

			const FBox aabb = box->get_aabb(nextSize);
			c += 6;

			//Z_XY
			dispatch_impl(res, mut, minc, c,
				std::rint(aabb.GetExtent().X),
				std::rint(aabb.GetExtent().Y),
				nextSize,
				std::rint(aabb.GetExtent().Z),
				EAxisPerm::Z_XY_0, b, co);

			//Z_YX
			dispatch_impl(res, mut, minc, c,
				std::rint(aabb.GetExtent().Y),
				std::rint(aabb.GetExtent().X),
				nextSize,
				std::rint(aabb.GetExtent().Z),
				EAxisPerm::Z_XY_1, b, co);

			//Y_XZ
			dispatch_impl(res, mut, minc, c,
				std::rint(aabb.GetExtent().X),
				std::rint(aabb.GetExtent().Z),
				nextSize,
				std::rint(aabb.GetExtent().Y),
				EAxisPerm::Y_XZ_0, b, co);

			//Y_ZX
			dispatch_impl(res, mut, minc, c,
				std::rint(aabb.GetExtent().Z),
				std::rint(aabb.GetExtent().X),
				nextSize,
				std::rint(aabb.GetExtent().Y),
				EAxisPerm::Y_XZ_1, b, co);

			//X_YZ
			dispatch_impl(res, mut, minc, c,
				std::rint(aabb.GetExtent().Y),
				std::rint(aabb.GetExtent().Z),
				nextSize,
				std::rint(aabb.GetExtent().X),
				EAxisPerm::X_YZ_0, b, co);

			//X_ZY
			dispatch_impl(res, mut, minc, c,
				std::rint(aabb.GetExtent().Z),
				std::rint(aabb.GetExtent().Y),
				nextSize,
				std::rint(aabb.GetExtent().X),
				EAxisPerm::X_YZ_1, b, co);

			while (c != 0) {}

			if (!isPacking_) break;

		}

		if (!isPacking_) break;

		if (res.empty()) {
			et++;
			mcc_++;
			if (useRandomBox && mcc_ > conf->MaxEmptryTries) break;
			if (useRandomBox && et > conf->EmptryTries) break;
			continue;
		}

		et = 0;

		std::sort(res.begin(), res.end(), [](const ::Detail::Result& _e1, const ::Detail::Result& _e2) {
			return _e1.weight < _e2.weight;
		});

		auto& r = res[0];

		const FBox tar = FBox(
			FVector(r.n0 - r.ext.X, r.n1 - r.ext.Y, r.h),
			FVector(r.n0 + r.ext.X, r.n1 + r.ext.Y, r.h + 2 * r.ext.Z));

		//DrawDebugBox(GetWorld(), tar.GetCenter(), tar.GetExtent(), FColor::Red, true);

		r.trans = ::Detail::make_transform(
			r.perm,
			tar,
			box->get_aabb(r.ext).GetCenter(),
			box->get_relative_location()
		);

		{
			std::lock_guard<std::mutex> lock(*m_);
			toSpawn_.push_back(r);
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
				const int32 i = n1 + n0 * N_;
				map_[i] = tar.Max.Z;
			}
		}

		tree_->recompute(mm);

	}

}//APacker::pack_impl

void APacker::Pack() {
	UWorld* world = GetWorld();
	if (!world) return;

	if (Config == nullptr) {
		std::cout << "Error: no Config set!" << std::endl;
		return;
	}

	Clear();

	using namespace MQT2;

	UPackerConfig* conf = Config->GetDefaultObject<UPackerConfig>();

	isPacking_ = true;
	std::cout << "start" << std::endl;
	std::cout << conf->Seed << std::endl;
	start_ = std::chrono::high_resolution_clock::now();
	future_ = std::make_unique<TFuture<bool>>(AsyncThread([this]() -> bool {
		pack_impl();
		const std::chrono::duration<double> ee = std::chrono::high_resolution_clock::now() - start_;
		std::cout << "done (" << ee.count() << "s)" << std::endl;
		last_time_ = ee.count();
		isPacking_ = false;
		return true;
		})
	);

}//APacker::Pack