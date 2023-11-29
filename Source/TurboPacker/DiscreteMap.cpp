
#include "DiscreteMap.h"

void ASpectralTester::Clear() {
	using namespace Util;

	UWorld* world = GetWorld();
	if (!world) return;

	FlushPersistentDebugLines(world);
}//ASpectralTester::Clear()

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

void ASpectralTester::TestHeightMap() {

	using namespace Util;

	UWorld* world = GetWorld();
	if (!world) return;

	FlushPersistentDebugLines(world);

	const FBox bounds(FVector(0.f), FVector(120, 80, 1700));
	DrawDebugBox(world, bounds.GetCenter(), bounds.GetExtent(), FColor::Blue, true);

	//1200*800*1700
	using Rollcage = HeightMap<120, 80, 1700>;
	Rollcage map;
	
	const FBox p(FVector(40, 20, 0), FVector(60, 60, 10));
	DrawDebugBox(world, p.GetCenter(), p.GetExtent(), FColor::Blue, true);

	map.push(p.GetCenter(), p.GetExtent());

	{
		const auto start = std::chrono::high_resolution_clock::now();
		const auto res = map.overlap(FVector(10.));
		const std::chrono::duration<double> ee = std::chrono::high_resolution_clock::now() - start;
		std::cout << "Overlap: " << ee.count() << "s" << std::endl;

		/*for(int32 y = 0; y < Rollcage::Y; ++y){
			for (int32 x = 0; x < Rollcage::X; ++x) {
				const int32 i = map.idx(x, y);
				if (i < 0 || i >= res[0].size()) continue;
				switch (res[0][i]) {
					case -1:
					DrawDebugPoint(world, FVector(x, y, map[{x, y}]), 1.f, FColor::Yellow, true);
					break;
					case 0:
					DrawDebugPoint(world, FVector(x, y, map[{x, y}]), 1.f, FColor::Green, true);
					break;
					case 1:
					DrawDebugPoint(world, FVector(x, y, map[{x, y}]), 1.f, FColor::Red, true);
					break;
				}
			}
		}*/

		for (int32 y = 0; y < Rollcage::Y; ++y) {
			for (int32 x = 0; x < Rollcage::X; ++x) {
				const int32 i = map.idx(x, y);
				if (i < 0 || i >= res[0].size()) continue;
				const int32 c = std::clamp(res[0][i] / Itensity, 0, 255);
				DrawDebugPoint(world, FVector(x, y, 0), 1.f, FColor(c, c, c), true);
			}
		}
	}

}//ASpectralTester::TestHeightMap