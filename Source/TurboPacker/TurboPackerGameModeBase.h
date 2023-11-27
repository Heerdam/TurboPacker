// Copyright Epic Games, Inc. All Rights Reserved.

#pragma once

#include "Util.h"

#include "GameFramework/GameModeBase.h"
#include "TurboPackerGameModeBase.generated.h"

/**
 * 
 */
UCLASS()
class TURBOPACKER_API ATurboPackerGameModeBase : public AGameModeBase {
	GENERATED_BODY()

public:
	ATurboPackerGameModeBase();
	void BeginPlay() override;
	
};//ATurboPackerGameModeBase
