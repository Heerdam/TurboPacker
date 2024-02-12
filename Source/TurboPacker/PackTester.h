#pragma once

#include "Util.h"
#include "GameFramework/Actor.h"
#include "GameFramework/Character.h"
#include "GameFramework/DefaultPawn.h"
#include "GameFramework/PlayerController.h"
#include "InputActionValue.h"

#include "MQT2.hpp"

#include "Async/Async.h"
#include "PackTester.generated.h"

/*
	Z_: Up axis
	_n: negativ up axis
	XY_: object axis on the world axis
	_N: rotation of N*90° in ccw around the up axis
*/
UENUM()
enum class EAxisPerm : uint8 {

	Z_XY_0, Z_XY_1, Z_XY_2, Z_XY_3,
	Y_XZ_0, Y_XZ_1, Y_XZ_2, Y_XZ_3,
	X_YZ_0, X_YZ_1, X_YZ_2, X_YZ_3,

	Z_n_XY_0, Z_n_XY_1, Z_n_XY_2, Z_n_XY_3,
	Y_n_XZ_0, Y_n_XZ_1, Y_n_XZ_2, Y_n_XZ_3,
	X_n_YZ_0, X_n_YZ_1, X_n_YZ_2, X_n_YZ_3

};//EAxisPerm

namespace Detail {

	inline FTransform make_transform(
		EAxisPerm _perm,
		const FBox& _target,
		const FVector& _pivot_offset,
		const FVector& _relative_offset
	);//make_transform

	struct Result {
		bool isRandomBox;
		FTransform trans;
		double weight;
		int32 n0, n1, h;
		int32 b_l, b_m, b_h;
		FVector ext;
		FVector ext_org;
		EAxisPerm perm;
		TSubclassOf<class APackerBox> box;
	};//Result

}//Detail

//-----------------------

/*
	ONLINE:		box is randomly picked/ generated. one after the other. 
				breaks when conditions is met or list is empty.
*/

UENUM(BlueprintType)
enum class PackMethod : uint8 {
	ONLINE
};//PackType

UENUM(BlueprintType)
enum class BoxGenerationType : uint8 {
	RANDOM, LIST
};//PackType

UENUM(BlueprintType)
enum class CostFunction : uint8 {
	SIMPLE, SIMPLE_HEIGHT
};//CostFunction

//-----------------------

USTRUCT(BlueprintType)
struct TURBOPACKER_API FPPair {
	GENERATED_BODY()

	UPROPERTY(EditAnywhere)
	FVector Size;

	UPROPERTY(EditAnywhere)
	int32 Count;
};//FPPair

//-----------------------

UCLASS(Blueprintable, BlueprintType)
class TURBOPACKER_API UBoxList : public UObject {
	GENERATED_BODY()
public:
	UPROPERTY(EditDefaultsOnly, Category = BoxList)
	TArray<FPPair> List;

	UPROPERTY(EditDefaultsOnly, Category = BoxList)
	bool ShuffleBoxes = false;
};//UBoxList

//-----------------------

UCLASS(Blueprintable, BlueprintType)
class TURBOPACKER_API UPackerConfig : public UObject {
	GENERATED_BODY()
public:

	UPROPERTY(EditDefaultsOnly, Category = General)
	bool MultiThreading = true;

	UPROPERTY(EditDefaultsOnly, Category = General)
	bool UseRandomSeed = true;

	UPROPERTY(EditDefaultsOnly, Category = General)
	uint64 Seed = 1234567890;

	UPROPERTY(EditDefaultsOnly, Category = General)
	FIntVector2 Bounds = FIntVector2(480);

	UPROPERTY(EditDefaultsOnly, Category = General)
	int32 Height = 150;

	UPROPERTY(EditDefaultsOnly, Category = General)
	PackMethod Method = PackMethod::ONLINE;

	UPROPERTY(EditDefaultsOnly, Category = General)
	BoxGenerationType BoxType = BoxGenerationType::RANDOM;

	UPROPERTY(EditDefaultsOnly, Category = General)
	CostFunction CostFunction = CostFunction::SIMPLE;

	//-------------------------

	UPROPERTY(EditDefaultsOnly, Category = Online)
	int32 EmptryTries = 25;

	UPROPERTY(EditDefaultsOnly, Category = Online)
	int32 MaxEmptryTries = 50;

	UPROPERTY(EditDefaultsOnly, Category = Online)
	int32 SetSize = 1;

	//-------------------------

	UPROPERTY(EditDefaultsOnly, Category = RandomBox)
	bool CubeRandomBoxes = true;

	UPROPERTY(EditDefaultsOnly, Category = RandomBox)
	double MinBoxVolume = 250.;

	UPROPERTY(EditDefaultsOnly, Category = RandomBox)
	double MaxBoxVolume = 1500.;

	//-------------------------

	UPROPERTY(EditDefaultsOnly, Category = BoxList)
	int32 ListIndex = 0;

	UPROPERTY(EditDefaultsOnly, Category = BoxList)
	TArray<TSubclassOf<UBoxList>> List;

};//UPackerConfig

UCLASS(Blueprintable, BlueprintType)
class TURBOPACKER_API APacker : public AActor {
	GENERATED_BODY()

	using Tree = MQT2::MedianQuadTree<int16, 15>;

	std::unique_ptr<std::mutex> m;
	std::vector< Detail::Result > toSpawn;

	std::queue<APackerBox*> q;

	bool isPacking = false;

	std::vector<int16> map;
	std::unique_ptr<Tree> tree;

	std::unique_ptr<TFuture<bool>> future;
	std::chrono::high_resolution_clock::time_point start;

	void pack_impl();

	[[nodiscard]] std::vector<::Detail::Result> overlap_impl(
		const int32 _ext0, 
		const int32 _ext1, 
		const int32 _h
	);
	
	void dispatch_impl(
		std::vector<::Detail::Result>& _res,
		std::mutex& _m,
		std::atomic<double>& _minc,
		std::atomic<int32>& _c,
		const int32 _n0, const int32 _n1,
		const int32 _h, EAxisPerm _perm,
		const TSubclassOf<APackerBox>& _b,
		const std::function<double(const Detail::Result&)>& _cost
	);

	double last_time = 0.;
	double vol = 0.;
	int32 bcc = 0;
	int32 mcc = 0;

public:

	UFUNCTION(BlueprintCallable)
	double get_pack_percent();

	UFUNCTION(BlueprintCallable)
	int32 get_bcc() { return bcc; }

	UFUNCTION(BlueprintCallable)
	int32 get_mcc() { return mcc; }

	UFUNCTION(BlueprintCallable)
	double get_time();
	
	//----------------

	UPROPERTY(EditAnywhere)
	bool AllowOverlap = true;

	UPROPERTY(EditDefaultsOnly)
	TSubclassOf<class UPackerConfig> Config;

	APacker();

	void Tick(float _delta) override;
	bool ShouldTickIfViewportsOnly() const override { return true; }

	UFUNCTION(BlueprintCallable, CallInEditor, Category = Packer)
	void Clear();

	UFUNCTION(BlueprintCallable, CallInEditor, Category = Packer)
	void Stop();

	UFUNCTION(BlueprintCallable, CallInEditor, Category = Packer)
	void Pack();


};//APacker

//--------------------------------------

UCLASS(Blueprintable, BlueprintType)
class TURBOPACKER_API AObserverController : public APlayerController {
	GENERATED_BODY()

	void Pack(const FInputActionValue& Value);
	void Clear(const FInputActionValue& Value);

public:
	UPROPERTY(EditDefaultsOnly)
	class UInputMappingContext* IA_context;

	UPROPERTY(EditDefaultsOnly)
	class APacker* Packer;

	AObserverController() = default;
	void BeginPlay() override;
	void OnPossess(APawn* aPawn) override;
	void SetupInputComponent() override;
};//AObserverController

//--------------------------------------

UCLASS(Blueprintable, BlueprintType)
class TURBOPACKER_API APackerBox : public AActor {
	GENERATED_BODY()

public:

	UPROPERTY(EditDefaultsOnly, Category = Boxler)
	UStaticMeshComponent* mesh = nullptr;

	APackerBox();
	virtual FBox get_aabb(const FVector _ext) const;
	FVector get_relative_location() const { return FVector(0.); }
	double get_weight() const;

};//APackerBox

UCLASS(Blueprintable, BlueprintType)
class TURBOPACKER_API ARandomBox : public APackerBox {
	GENERATED_BODY()

public:

	UPROPERTY(EditDefaultsOnly, Category = Boxler)
	double MaxSize = 1500.;

	ARandomBox();

	virtual FBox get_aabb(const FVector _ext) const override;
	void set_to_size(const FVector& _new_size, const double _min_size, const double _max_size);

};//ARandomBox
