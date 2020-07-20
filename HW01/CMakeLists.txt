cmake_minimum_required(VERSION 3.16)
project(HW01)

set(CMAKE_CXX_STANDARD 11)

aux_source_directory(. DIR_SRC)

add_executable(${PROJECT_NAME} ${DIR_SRC})

target_link_libraries(${PROJECT_NAME} ${LIB_SURFACE_MESH} ${OPENMP})