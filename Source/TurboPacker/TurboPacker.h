// Copyright Epic Games, Inc. All Rights Reserved.

#pragma once

#include "CoreMinimal.h"

class FTurboPackerGameModule : public IModuleInterface {
	void* handle;
public:
	void StartupModule() override;
	void ShutdownModule() override;
	bool IsGameModule() const override { return true; }
};//FTurboPackerGameModule
