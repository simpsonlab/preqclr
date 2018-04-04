cmake_minimum_required(VERSION 3.2)
project(spoa)

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/bin)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -pedantic -march=native")
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

option(spoa_build_executable "Build spoa standalone tool" OFF)
option(spoa_build_tests "Build spoa unit tests" OFF)

add_library(spoa STATIC
    src/alignment_engine.cpp
    src/graph.cpp
    src/simd_alignment_engine.cpp
    src/sisd_alignment_engine.cpp)

target_include_directories(spoa PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>)

install(TARGETS spoa DESTINATION lib)
install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/include/spoa DESTINATION include)

if (spoa_build_executable)
    add_executable(spoa_bin
        src/sequence.cpp
        src/main.cpp)

    if (NOT TARGET bioparser)
        add_subdirectory(vendor/bioparser EXCLUDE_FROM_ALL)
    endif()

    target_link_libraries(spoa_bin spoa bioparser)
    set_target_properties(spoa_bin PROPERTIES OUTPUT_NAME spoa)

    install(TARGETS spoa_bin DESTINATION bin)
endif(spoa_build_executable)

if (spoa_build_tests)
    set(spoa_test_data_path ${PROJECT_SOURCE_DIR}/test/data/)
    configure_file(${PROJECT_SOURCE_DIR}/test/spoa_test_config.h.in
        ${PROJECT_BINARY_DIR}/config/spoa_test_config.h)
    include_directories(${PROJECT_BINARY_DIR}/config)
    include_directories(${PROJECT_SOURCE_DIR}/src)

    add_executable(spoa_test
        src/sequence.cpp
        test/spoa_test.cpp)

    if (NOT TARGET bioparser)
        add_subdirectory(vendor/bioparser EXCLUDE_FROM_ALL)
    endif()
    add_subdirectory(vendor/googletest/googletest EXCLUDE_FROM_ALL)

    target_link_libraries(spoa_test spoa bioparser gtest_main)
endif(spoa_build_tests)
