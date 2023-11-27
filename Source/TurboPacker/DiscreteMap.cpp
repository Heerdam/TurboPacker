
#include "DiscreteMap.h"

void ASpectralTester::TestInfGrid() {

	using namespace Util;

	UWorld* world = GetWorld();
	if (!world) return;

	FlushPersistentDebugLines(world);

	std::random_device rd;
	Rand rand(rd());

	Eigen::Vector3d size(10., 10., 10.);

	InfBucketGrid<bool, FVector> map(size);

	Dist d(0., 50.);
	for (int32 i = 0; i < 100000; ++i) {
		map.push(true, FVector(d(rand), d(rand), d(rand)));
	}

	for (int32 z = 0; z < 5; ++z) {
		for (int32 y = 0; y < 5; ++y) {
			for (int32 x = 0; x < 5; ++x) {

				DrawDebugBox(world, FVector(x * 10.f, y * 10.f, z * 10.f), FVector(4.9f), FColor::Blue, true);

			}
		}
	}

	auto b1 = map.radial_search(FVector(25.f), 5.f);

	std::cout << b1.size() << std::endl;

	for (const auto t : b1) {
		const auto& [tt, p] = *t;
		DrawDebugPoint(world, p, 1.f, FColor::Red, true);
	}

}//ASpectralTester::TestInfGrid