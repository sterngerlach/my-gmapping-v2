
# CMakeLists.txt

cmake_minimum_required(VERSION 3.8)

project(my_gmapping CXX)

# Find libpng
find_package(PNG REQUIRED)

# Find Boost package
set (BOOST_USE_STATIC_LIBS ON)
find_package(Boost 1.65.1 REQUIRED system timer)

# Find Eigen package
find_package(Eigen3 3.3 REQUIRED NO_MODULE)

# Find Neon intrinsics
list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake/Modules)
find_package(ARM)

# Find OpenMP
find_package(OpenMP REQUIRED)

set(CMAKE_DEBUG_POSTFIX "d")

add_library(gmapping)

add_subdirectory(include)
add_subdirectory(src)

target_include_directories(gmapping PUBLIC ${PROJECT_SOURCE_DIR}/include)
target_include_directories(gmapping PUBLIC ${EIGEN3_INCLUDE_DIR})
target_include_directories(gmapping PUBLIC ${Boost_INCLUDE_DIRS})
target_link_libraries(gmapping PUBLIC m rt pthread dl)
target_link_libraries(gmapping PUBLIC Eigen3::Eigen)
target_link_libraries(gmapping PUBLIC PNG::PNG)
target_link_libraries(gmapping PUBLIC ${Boost_LIBRARIES})

if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    if (CMAKE_CXX_COMPILER_VERSION VERSION_LESS 9.1)
        target_link_libraries(gmapping PUBLIC stdc++fs)
    endif()
endif()

target_compile_options(gmapping PUBLIC
    $<$<CXX_COMPILER_ID:GNU>:-pedantic -Wall -Wextra>
    $<$<CONFIG:Release>:-O3>
    $<$<AND:$<CXX_COMPILER_ID:GNU>,$<CONFIG:Debug>>:-O0 -ggdb>)
target_link_libraries(gmapping PUBLIC OpenMP::OpenMP_CXX)
target_compile_features(gmapping PUBLIC cxx_std_17)
set_target_properties(gmapping PROPERTIES
    DEBUG_POSTFIX ${CMAKE_DEBUG_POSTFIX})

if (NEON_TRUE)
    message(STATUS "Neon found with compiler flag: -mfpu=neon -D__NEON__")
    target_compile_options(gmapping PUBLIC -mfpu=neon -D__NEON__)
    target_compile_options(gmapping PUBLIC -munaligned-access)
endif (NEON_TRUE)

add_executable(slam_launch
    ${PROJECT_SOURCE_DIR}/src/my_gmapping/slam_launcher.cpp)
target_link_libraries(slam_launch PUBLIC gmapping)
set_target_properties(slam_launch PROPERTIES
    DEBUG_POSTFIX ${CMAKE_DEBUG_POSTFIX})
