
# CMakeLists.txt

add_subdirectory(carmen)

target_sources(gmapping PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/gnuplot_helper.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/map_saver.cpp)
