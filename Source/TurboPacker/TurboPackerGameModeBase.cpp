
#include "TurboPackerGameModeBase.h"

ATurboPackerGameModeBase::ATurboPackerGameModeBase() {
	PrimaryActorTick.bCanEverTick = true;
	LStream::set();
}

void ATurboPackerGameModeBase::BeginPlay() {

	Eigen::FFT<float> fft;

	std::vector<float> timevec;
	for (int32 i = 0; i < 1000; ++i)
		timevec.push_back(float(i));
	std::vector<std::complex<float>> freqvec;

	fft.fwd(freqvec, timevec);
	fft.inv(timevec, freqvec);

	//UE_LOG(LogTemp, Warning, TEXT("Ding"));

}//ATurboPackerGameModeBase::BeginPlay