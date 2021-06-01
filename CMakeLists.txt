cmake_minimum_required(VERSION 3.13)
project(gdon3d VERSION 1.0.0 LANGUAGES C CXX)

set(CMAKE_C_STANDARD 99)
set(CMAKE_CXX_STANDARD 11)

# Out-of-source builds guard
file(TO_CMAKE_PATH "${PROJECT_BINARY_DIR}/CMakeLists.txt" LOC_PATH)
if(EXISTS "${LOC_PATH}")
    message(FATAL_ERROR "You cannot build in a source directory (or any directory with a CMakeLists.txt file). Please make a build subdirectory. Feel free to remove CMakeCache.txt and CMakeFiles.")
endif()

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Choose Release or Debug" FORCE)
endif()

set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/tools/cmake/modules)

# External dep
find_package(OpenMP REQUIRED)
find_package(SuiteSparse REQUIRED)
find_package(HDF5 REQUIRED)

add_library(gdon3d STATIC "")

target_compile_definitions(gdon3d PUBLIC _DOUBLE_PRECISION)

target_include_directories(gdon3d PUBLIC ${SuiteSparse_INCLUDE_DIRS}
                                         ${HDF5_INCLUDE_DIRS})

target_link_libraries(gdon3d PUBLIC ${SuiteSparse_LIBRARIES}  
                                    ${HDF5_C_LIBRARIES}
                                    OpenMP::OpenMP_C
                                    m)

target_compile_options(gdon3d PRIVATE "${OpenMP_C_FLAGS}" -g)
add_subdirectory(src)
add_subdirectory(examples)