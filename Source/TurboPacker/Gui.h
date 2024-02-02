#pragma once

#include "Util.h"

#include "EditorUtilityWidget.h"

#include "Gui.generated.h"

UCLASS(Blueprintable, BlueprintType)
class TURBOPACKER_API UOnlinePackerWidget : public UEditorUtilityWidget{
	GENERATED_BODY()

public:

	UPROPERTY(EditDefaultsOnly)
	class AOnlinePacker* Packer;

	UFUNCTION()
	double GetDensity();

	UFUNCTION()
	int32 GetPackedBoxes();

	UFUNCTION()
	int32 GetMissedBoxes();
	
};//UOnlinePackerWidget
