#pragma once

#include "CoreMinimal.h"

#include <sstream>
#include <iostream>

class TURBOPACKER_API LStream : public std::stringbuf {

	static LStream* ptr;

	LStream() {
		set();
	}

protected:
	int sync() {
		UE_LOG(LogTemp, Warning, TEXT("%s"), *FString(str().c_str()));
		str("");
		return std::stringbuf::sync();
	}

public:
	static void set() {
		std::cout.rdbuf(ptr);
	}
};//LStream
