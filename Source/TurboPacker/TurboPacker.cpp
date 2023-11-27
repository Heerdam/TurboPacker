// Copyright Epic Games, Inc. All Rights Reserved.

#include "TurboPacker.h"
#include "Modules/ModuleManager.h"

void FTurboPackerGameModule::StartupModule() {
	{
		const auto path = FPaths::ProjectDir().Append("Libs/FFTW/libfftw3f-3.dll");
		handle1 = FPlatformProcess::GetDllHandle(*path);
	}
	{
		const auto path = FPaths::ProjectDir().Append("Libs/FFTW/libfftw3l-3.dll");
		handle2 = FPlatformProcess::GetDllHandle(*path);
	}
	{
		const auto path = FPaths::ProjectDir().Append("Libs/FFTW/libfftw3-3.dll");
		handle3 = FPlatformProcess::GetDllHandle(*path);
	}
	//UE_LOG(LogTemp, Warning, TEXT("Project Directory: %s"), *path);
}//FTurboPackerGameModule::StartupModule

void FTurboPackerGameModule::ShutdownModule() {
	FPlatformProcess::FreeDllHandle(handle1);
	FPlatformProcess::FreeDllHandle(handle2);
	FPlatformProcess::FreeDllHandle(handle3);
}//FTurboPackerGameModule::ShutdownModule

IMPLEMENT_PRIMARY_GAME_MODULE(FTurboPackerGameModule, TurboPacker, "TurboPacker" );
