#pragma once

#include "Util.h"
#include "GameFramework/Actor.h"

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
	);

}//Detail

UCLASS(Blueprintable, BlueprintType)
class TURBOPACKER_API APackTester : public AActor {
	GENERATED_BODY()

	int32 idx = 0;
	TArray<TPair<FTransform, int32>> cache;

public:

	UPROPERTY(EditAnywhere)
	FIntVector Bounds;

	UPROPERTY(EditDefaultsOnly)
	TArray<TSubclassOf<class APackerBox>> Boxes;

	UFUNCTION(BlueprintCallable, CallInEditor, Category = Packer)
	void Clear();

	UFUNCTION(BlueprintCallable, CallInEditor, Category = Packer)
	void Pack();

	UFUNCTION(BlueprintCallable, CallInEditor, Category = Packer)
	void StepForward();

	UFUNCTION(BlueprintCallable, CallInEditor, Category = Packer)
	void StepBack();

};//ASpectralPacker

UCLASS(Blueprintable, BlueprintType)
class TURBOPACKER_API APackerBox : public AActor {
	GENERATED_BODY()

public:

	UPROPERTY(EditDefaultsOnly, Category = Boxler)
	UStaticMeshComponent* mesh = nullptr;

	APackerBox();
	FBox get_aabb() const;
	virtual FVector get_relative_location() const { return FVector(0.); }

};//APackerBox