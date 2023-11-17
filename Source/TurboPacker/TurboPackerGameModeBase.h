// Copyright Epic Games, Inc. All Rights Reserved.

#pragma once

#include <Eigen/Core>
#include <unsupported/Eigen/FFT>


#include "CoreMinimal.h"
#include "GameFramework/GameModeBase.h"
#include "TurboPackerGameModeBase.generated.h"

/**
 * 
 */
UCLASS()
class TURBOPACKER_API ATurboPackerGameModeBase : public AGameModeBase {
	GENERATED_BODY()

public:

	void BeginPlay() override;
	
};//ATurboPackerGameModeBase
