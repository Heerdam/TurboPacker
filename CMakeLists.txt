
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

option(BUILD_AS_LIBRARY "Use as library" FALSE)

if (NOT BUILD_AS_LIBRARY)

	add_executable( TP_App src/main.cpp )

	# --------------- RAYLIB ---------------
	message(STATUS "Fetching Raylib...")
	FetchContent_Declare(
		RAYLIB
		GIT_REPOSITORY  https://github.com/raysan5/raylib.git
		GIT_TAG         5.0
	)
	option(BUILD_EXAMPLES "Build example of how to use this library." OFF)
	option(BUILD_SHARED_LIBS "Build shared libs." OFF)
	set(GRAPHICS "GRAPHICS_API_OPENGL_43" FORCE)
	FetchContent_MakeAvailable(RAYLIB)

	# ------------------ IMGUI ------------------
	message(STATUS "Fetching imgui...")
	FetchContent_Declare(
		imgui
		GIT_REPOSITORY      https://github.com/ocornut/imgui.git
		GIT_TAG             v1.90.4
	)
	FetchContent_GetProperties(imgui)
	if(NOT imgui_POPULATED)
		FetchContent_Populate(imgui)
		set(IMGUI_FILES
			"${imgui_SOURCE_DIR}/imgui_demo.cpp"
			"${imgui_SOURCE_DIR}/imgui_draw.cpp"
			"${imgui_SOURCE_DIR}/imgui_widgets.cpp"
			"${imgui_SOURCE_DIR}/imgui_tables.cpp"
			"${imgui_SOURCE_DIR}/imgui.cpp"		
		)
		add_library(imgui "${IMGUI_FILES}")
		target_include_directories(TP_App PUBLIC "${imgui_SOURCE_DIR}")
	endif()

	# ------------------ RAYGUI ------------------
	message(STATUS "Fetching raygui...")
	FetchContent_Declare(
		rimgui
		GIT_REPOSITORY      https://github.com/raylib-extras/rlImGui.git
		GIT_TAG 			main
	)
	FetchContent_GetProperties(rimgui)
	if(NOT Rimgui_POPULATED)
		FetchContent_Populate(rimgui)
		set(RIMGUI_FILES
			"${rimgui_SOURCE_DIR}/rlImGui.cpp"
		)
		add_library(rimgui "${RIMGUI_FILES}")
		target_include_directories(TP_App PUBLIC "${rimgui_SOURCE_DIR}")
	endif()
	target_include_directories(rimgui PUBLIC "${imgui_SOURCE_DIR}")
	target_link_libraries(rimgui PUBLIC imgui raylib)

	# ------------------ LODEPNG ------------------
	message(STATUS "Fetching lodepng...")
	FetchContent_Declare(
		lodepng
		GIT_REPOSITORY      https://github.com/lvandeve/lodepng.git
	)
	FetchContent_GetProperties(lodepng)
	if(NOT lodepng_POPULATED)
		FetchContent_Populate(lodepng)

		set(LODEPNG_FILES
			"${lodepng_SOURCE_DIR}/lodepng.cpp"	
		)

		add_library(lodepng "${LODEPNG_FILES}")
		target_include_directories(TP_App PUBLIC "${lodepng_SOURCE_DIR}")
	endif()

	# ------------------ JSON ------------------
	message(STATUS "Fetching json...")
	FetchContent_Declare(json 
		URL https://github.com/nlohmann/json/releases/download/v3.11.3/json.tar.xz
		DOWNLOAD_EXTRACT_TIMESTAMP TRUE
	)
	FetchContent_MakeAvailable(json)

endif()

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
	GIT_TAG             1.0.1
)
FetchContent_MakeAvailable(glm)

if (NOT BUILD_AS_LIBRARY)
	target_include_directories(TP_App PUBLIC include )
	target_link_libraries(TP_App MQT2 glm::glm raylib imgui rimgui lodepng nlohmann_json::nlohmann_json)
endif()

# -----------------------------------------

target_include_directories(TP INTERFACE include )
target_link_libraries(TP INTERFACE MQT2 glm::glm )

