cmake_minimum_required(VERSION 3.16)
project(DGP_Homework)

set(CMAKE_CXX_STANDARD 11)

# added -fopenmp
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Xpreprocessor -fopenmp")

set(LIB_SURFACE_MESH /usr/local/lib/libSurfaceMesh.1.0.dylib)
set(OPENMP /usr/local/lib/libomp.dylib)
set(INCLUDE /usr/local/include)

include_directories(${INCLUDE})

find_package(Eigen3 3.3.0)
if (Eigen3_FOUND)
    include_directories(${EIGEN3_INCLUDE_DIR})
else ()
    message(can't find eigen3)
    return()
endif ()

aux_source_directory(. DIR_SRC)
aux_source_directory(./utils UTILS_SRC)

add_subdirectory(HW01)

if (MSVC)
    target_compile_definitions(${PROJECT_NAME} PRIVATE _USE_MATH_DEFINES)
endif (MSVC)
