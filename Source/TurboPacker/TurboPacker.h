// Copyright Epic Games, Inc. All Rights Reserved.

#pragma once

#include "CoreMinimal.h"

class FTurboPackerGameModule : public IModuleInterface {
	void* handle1;
	void* handle2;
	void* handle3;
public:
	void StartupModule() override;
	void ShutdownModule() override;
	bool IsGameModule() const override { return true; }
};//FTurboPackerGameModule
