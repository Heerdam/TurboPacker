// Copyright Epic Games, Inc. All Rights Reserved.

using System.IO;
using UnrealBuildTool;

public class TurboPacker : ModuleRules
{
	public TurboPacker(ReadOnlyTargetRules Target) : base(Target)
	{
		PCHUsage = PCHUsageMode.UseExplicitOrSharedPCHs;
	
		PublicDependencyModuleNames.AddRange(new string[] { 
            "Core", 
            "CoreUObject", 
            "Engine", 
            "InputCore",
            "EnhancedInput",
            "OutputLog" 
        });

		PrivateDependencyModuleNames.AddRange(new string[] {  });

        PrivateIncludePaths.AddRange(
            new string[] {
                Path.Combine(ModuleDirectory, "../../Ext/robin/include"),
                Path.Combine(ModuleDirectory, "../../Ext/libmorton/include"),
                Path.Combine(ModuleDirectory, "../../Libs/eigen"),
                Path.Combine(ModuleDirectory, "../../Libs/FFTW"),
                Path.Combine(ModuleDirectory, "../../Libs/mqt/include")

            }
        );

        PublicDefinitions.Add("__LDBL_MANT_DIG__=53");
        PublicDefinitions.Add("EIGEN_FFTW_DEFAULT");
       // PublicDefinitions.Add("__x86_64__");
        //PublicDefinitions.Add("__GNUC__=3");
        //PublicDefinitions.Add("__GNUC_MINOR__=6");

        PublicDelayLoadDLLs.Add("libfftw3-3.dll");
        PublicDelayLoadDLLs.Add("libfftw3f-3.dll");
        PublicDelayLoadDLLs.Add("libfftw3l-3.dll");

        PublicAdditionalLibraries.Add(Path.Combine(ModuleDirectory, "../../Libs/FFTW/libfftw3-3.lib"));
        PublicAdditionalLibraries.Add(Path.Combine(ModuleDirectory, "../../Libs/FFTW/libfftw3f-3.lib"));
        PublicAdditionalLibraries.Add(Path.Combine(ModuleDirectory, "../../Libs/FFTW/libfftw3l-3.lib"));

        //AddEngineThirdPartyPrivateStaticDependencies(Target, "Eigen");

        // Uncomment if you are using Slate UI
        // PrivateDependencyModuleNames.AddRange(new string[] { "Slate", "SlateCore" });

        // Uncomment if you are using online features
        // PrivateDependencyModuleNames.Add("OnlineSubsystem");

        // To include OnlineSubsystemSteam, add it to the plugins section in your uproject file with the Enabled attribute set to true
    }
}
