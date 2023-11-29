
#include "TurboPackerGameModeBase.h"

ATurboPackerGameModeBase::ATurboPackerGameModeBase() {
	PrimaryActorTick.bCanEverTick = true;
	LStream::set();
}

void ATurboPackerGameModeBase::BeginPlay() {

}//ATurboPackerGameModeBase::BeginPlay