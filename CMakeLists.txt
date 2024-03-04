
cmake_minimum_required(VERSION 3.16)

project(TP)

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

include(FetchContent)

add_library( TP INTERFACE )

add_executable(TPtest
   src/main.cpp
)

# ------------------ MQT ------------------
message(STATUS "Fetching mqt...")
FetchContent_Declare(
	MQT2
	GIT_REPOSITORY https://github.com/Heerdam/MedianQuadtree.git
)
FetchContent_MakeAvailable(MQT2)

# ------------------ GLM ------------------
message(STATUS "Fetching glm...")
FetchContent_Declare(
    glm
    GIT_REPOSITORY      https://github.com/g-truc/glm.git
    GIT_TAG             1.0.0
)
FetchContent_MakeAvailable(glm)

# -----------------------------------------

target_include_directories(TP INTERFACE include )
target_include_directories(TPtest PUBLIC include )

target_link_libraries(TP INTERFACE MQT2 glm::glm )
target_link_libraries(TPtest MQT2 glm::glm )
