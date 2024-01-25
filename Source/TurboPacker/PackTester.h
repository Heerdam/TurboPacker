#pragma once

#include "Util.h"
#include "GameFramework/Actor.h"
#include "lode/lodepng.h"

#include "MQT2.hpp"

#include "Async/Async.h"
#include "PackTester.generated.h"

/*
	Z_: Up axis
	_n: negativ up axis
	XY_: object axis on the world axis
	_N: rotation of N*90� in ccw around the up axis
*/
UENUM()
enum class EAxisPerm : uint8 {

	Z_XY_0, Z_XY_1, Z_XY_2, Z_XY_3,
	Y_XZ_0, Y_XZ_1, Y_XZ_2, Y_XZ_3,
	X_YZ_0, X_YZ_1, X_YZ_2, X_YZ_3,

	Z_n_XY_0, Z_n_XY_1, Z_n_XY_2, Z_n_XY_3,
	Y_n_XZ_0, Y_n_XZ_1, Y_n_XZ_2, Y_n_XZ_3,
	X_n_YZ_0, X_n_YZ_1, X_n_YZ_2, X_n_YZ_3

};//EAxisPerm

namespace Detail {

	inline FTransform make_transform(
		EAxisPerm _perm,
		const FBox& _target,
		const FVector& _pivot_offset,
		const FVector& _relative_offset
	);//make_transform

	struct Result {
		double weight;
		int32 n0, n1, h;
		FVector ext;
		EAxisPerm perm;
		TSubclassOf<class APackerBox> box;
	};//Result

}//Detail

//-----------------------

UCLASS(Blueprintable, BlueprintType)
class TURBOPACKER_API APackTester : public AActor {
	GENERATED_BODY()

	using Tree = MQT2::MedianQuadTree<float, 15>;

	std::unique_ptr<std::mutex> m;
	std::vector<std::pair<TSubclassOf<APackerBox>, FTransform>> toSpawn;

	bool isPacking = false;

	std::vector<float> map;
	std::unique_ptr<Tree> tree;

	std::unique_ptr<TFuture<bool>> future;

	void pack_impl();

public:

	UPROPERTY(EditAnywhere)
	int32 Bounds = 480;

	UPROPERTY(EditAnywhere)
	int32 Height = 150;

	UPROPERTY(EditDefaultsOnly)
	TArray<TSubclassOf<class APackerBox>> Boxes;

	APackTester();

	void Tick(float _delta) override;
	virtual bool ShouldTickIfViewportsOnly() const override { return true; }

	UFUNCTION(BlueprintCallable, CallInEditor, Category = Packer)
	void Clear();

	UFUNCTION(BlueprintCallable, CallInEditor, Category = Packer)
	void Pack();

};//APackTester

//-----------------------
USTRUCT(BlueprintType)
struct TURBOPACKER_API FPPair {
	GENERATED_BODY()

	UPROPERTY(EditAnywhere)
	TSubclassOf<class APackerBox> Type;

	UPROPERTY(EditAnywhere)
	int32 Count;
};//FPPair

UCLASS(Blueprintable, BlueprintType)
class TURBOPACKER_API UOnlinePackerConfig : public UObject {
	GENERATED_BODY()
public:
	UPROPERTY(EditDefaultsOnly)
	int32 Bounds = 480;

	UPROPERTY(EditDefaultsOnly)
	int32 Height = 150;

	UPROPERTY(EditDefaultsOnly)
	TArray<FPPair> Boxes;
};//UOnlinePackerConfig

UCLASS(Blueprintable, BlueprintType)
class TURBOPACKER_API AOnlinePacker : public AActor {
	GENERATED_BODY()

	using Tree = MQT2::MedianQuadTree<float, 15>;

	std::unique_ptr<std::mutex> m;
	std::vector< std::pair<TSubclassOf<APackerBox>, FTransform> > toSpawn;

	std::queue<APackerBox*> q;

	bool isPacking = false;

	std::vector<float> map;
	std::unique_ptr<Tree> tree;

	std::unique_ptr<TFuture<bool>> future;
	std::chrono::high_resolution_clock::time_point start;

	void pack_impl();

	double vol = 0.;

public:

	UPROPERTY(EditDefaultsOnly)
	TSubclassOf<class UOnlinePackerConfig> Config;

	AOnlinePacker();

	void Tick(float _delta) override;
	virtual bool ShouldTickIfViewportsOnly() const override { return true; }

	UFUNCTION(BlueprintCallable, CallInEditor, Category = Packer)
	void Clear();

	UFUNCTION(BlueprintCallable, CallInEditor, Category = Packer)
	void Pack();

};//AOnlinePacker

//--------------------------------------

UCLASS(Blueprintable, BlueprintType)
class TURBOPACKER_API APackerBox : public AActor {
	GENERATED_BODY()

public:

	UPROPERTY(EditDefaultsOnly, Category = Boxler)
	UStaticMeshComponent* mesh = nullptr;

	APackerBox();
	FBox get_aabb() const;
	virtual FVector get_relative_location() const { return FVector(0.); }
	double get_weight() const;

};//APackerBox

//--------------------------------------

template<class R, bool LogScaling>
inline void image_real(
	const int32 _n0,
	const int32 _n1,
	R* _buffer,
	const std::string& _filename,
	const std::string& _folder = "output"
) {
	{
		const std::filesystem::path path = std::filesystem::path(TCHAR_TO_UTF8(*FPaths::ProjectDir())) / _folder;
		if (!std::filesystem::exists(path))
			std::filesystem::create_directory(path);
	}

	double max = -std::numeric_limits<double>::infinity();
	for (int32 i = 0; i < _n0 * _n1; ++i)
		max = std::max<double>(max, std::abs(_buffer[i]));

	std::vector<unsigned char> img;
	for (int32 n0 = 0; n0 < _n0; ++n0) {
		for (int32 n1 = 0; n1 < _n1; ++n1) {
			const int32 i1 = n1 + n0 * _n1;
			unsigned char c = 0;

			if constexpr (LogScaling) {
				c = std::abs(1. - std::abs(std::log(std::abs(1. + _buffer[i1]))) / std::log(max));
			} else {
				const double frac = 1. / max;
				c = 255 * (std::abs(_buffer[i1]) * frac);
			}

			img.push_back(c);
			img.push_back(c);
			img.push_back(c);
			img.push_back(255);
		}
	}
	const std::filesystem::path path = std::filesystem::path(TCHAR_TO_UTF8(*FPaths::ProjectDir())) / _folder / _filename;
	lodepng::encode(path.string(), img.data(), _n1, _n0);
}//image_real