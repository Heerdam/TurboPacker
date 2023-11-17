// Copyright Epic Games, Inc. All Rights Reserved.

#include "TurboPacker.h"
#include "Modules/ModuleManager.h"

void FTurboPackerGameModule::StartupModule() {
	const auto path = FPaths::ProjectDir().Append("Libs/FFTW/libfftw3f-3.dll");
	handle = FPlatformProcess::GetDllHandle(*path);
	//UE_LOG(LogTemp, Warning, TEXT("Project Directory: %s"), *path);
}//FTurboPackerGameModule::StartupModule

void FTurboPackerGameModule::ShutdownModule() {
	FPlatformProcess::FreeDllHandle(handle);
}//FTurboPackerGameModule::ShutdownModule

IMPLEMENT_PRIMARY_GAME_MODULE(FTurboPackerGameModule, TurboPacker, "TurboPacker" );
