
cmake_minimum_required(VERSION 3.5)

project(MQT2)
project(MQT2 LANGUAGES CUDA CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
if (MSVC)
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")
	set(CMAKE_CXX_FLAGS_DEBUG)
	set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} /O2 /DNDEBUG")
else()
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")
	set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -g -ggdb -fsanitize=address,undefined -pedantic" )
	set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3 -march=native -DNDEBUG -funroll-loops -mavx2")
endif()

add_library( MQT2 INTERFACE )

add_executable(mqttest
   src/main.cpp
   src/MQT2.cu

)

#-ffast-math

set_target_properties(mqttest PROPERTIES CUDA_ARCHITECTURES "61")

target_include_directories(MQT2 INTERFACE include )
target_include_directories(mqttest PUBLIC include )

target_link_libraries(MQT2 INTERFACE)
target_link_libraries(mqttest)
