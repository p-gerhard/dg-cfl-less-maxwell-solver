cmake_minimum_required(VERSION 3.13)

target_sources(gdon3d PUBLIC gdon3d.h)

target_include_directories(gdon3d PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

add_subdirectory(dg)
add_subdirectory(io)
add_subdirectory(maths)
add_subdirectory(mesh)
add_subdirectory(models)
add_subdirectory(utils)
add_subdirectory(simulation)