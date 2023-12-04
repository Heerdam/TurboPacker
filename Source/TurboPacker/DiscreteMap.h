#pragma once

#include "Util.h"

#include "GameFramework/Actor.h"

#include "DiscreteMap.generated.h"

//1200*800*1700 - 1,632,000,000
//  - 623,208,300

namespace TurboPacker {

	namespace Spectral {

		namespace Detail {

		}// Detail

		class TURBOPACKER_API SpectralMap {

			Util::InfBucketGrid<bool, FVector> grid;

		public:



		}; //SpectralMap

	}//Spectral

}//TurboPacker

UCLASS(Blueprintable, BlueprintType)
class TURBOPACKER_API ASpectralTester : public AActor {
	GENERATED_BODY()

public:
	UFUNCTION(BlueprintCallable, CallInEditor, Category = Test)
	void Clear();

	UFUNCTION(BlueprintCallable, CallInEditor, Category = Test)
	void TestInfGrid();

	UFUNCTION(BlueprintCallable, CallInEditor, Category = Test)
	void TestHeightMap();

	UFUNCTION(BlueprintCallable, CallInEditor, Category = Test)
	void TestFFTW();

	UPROPERTY(EditAnywhere)
	int32 Itensity = 8;

};