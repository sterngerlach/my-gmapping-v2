
# CMakeLists.txt

target_sources(gmapping PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/bresenham.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/memory_usage.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/util.cpp)

add_subdirectory(grid_map_new)
add_subdirectory(hw)
add_subdirectory(io)
add_subdirectory(mapping)
add_subdirectory(metric)
