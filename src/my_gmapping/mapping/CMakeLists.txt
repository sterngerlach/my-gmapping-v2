
# CMakeLists.txt

target_sources(gmapping PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/covariance_estimator.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/grid_map_builder.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/likelihood_greedy_endpoint.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/map_builder.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/motion_model_relative_pose.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/particle.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/scan_interpolator.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/scan_matcher_correlative.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/scan_matcher_correlative_fpga.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/scan_matcher_hill_climbing.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/weight_normalizer.cpp)
