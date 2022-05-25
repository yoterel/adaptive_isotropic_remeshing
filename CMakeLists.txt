cmake_minimum_required(VERSION 3.10)

project(adaptive_remesh VERSION 0.1.0)
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED True)
set(CMAKE_BUILD_TYPE "Release")  # change to Debug for debugging
#set(CMAKE_POSITION_INDEPENDENT_CODE ON)
set(EXT_DIR ${CMAKE_CURRENT_SOURCE_DIR}/ext)
#set(NPE_PYTHON_EXECUTABLE "/path/to/virtual/env/or/any/python/bin")  # uncomment to set a path to a specific python release (default is system python version)

find_package(Threads REQUIRED)
find_package(Eigen3 REQUIRED)
configure_file(version.h.in version.h)

include_directories(
	# libigl
	${EXT_DIR}/libigl/include
	)

list(APPEND CMAKE_MODULE_PATH ${EXT_DIR}/numpyeigen/cmake)
include(numpyeigen)
	
add_library(adaptive_remesh_lib
	SHARED

	# Headers
	src/collapse_edges.h
	src/equalize_valences.h
	src/adaptive_remesh_botsch.h
	src/split_edges.h
	src/split_edges_until_bound.h
	src/tangential_smoothing.h

	# Source
	src/collapse_edges.cpp
	src/equalize_valences.cpp
	src/adaptive_remesh_botsch.cpp
	src/split_edges.cpp
	src/split_edges_until_bound.cpp
	src/tangential_smoothing.cpp

	)

target_link_libraries(adaptive_remesh_lib Threads::Threads Eigen3::Eigen)

add_executable(${PROJECT_NAME} main.cpp)

target_link_libraries(${PROJECT_NAME} adaptive_remesh_lib)

target_include_directories(${PROJECT_NAME} PUBLIC
                           ${PROJECT_BINARY_DIR}
                           )

npe_add_module(py_ada_remesh
				BINDING_SOURCES
				src/bindings.cpp
				)

target_include_directories(py_ada_remesh PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/src)
target_link_libraries(py_ada_remesh PUBLIC adaptive_remesh_lib)
