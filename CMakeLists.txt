cmake_minimum_required(VERSION 3.2)
project(preqclrhtzy C CXX)

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/bin)

#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -pedantic")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -O3")
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED OFF)
set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS}")

find_package(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

option(racon_build_tests "Build preqclrhtzy unit tests" OFF)
option(racon_build_wrapper "Build preqclrhtzy wrapper" OFF)

add_executable(preqclrhtzy
    #src/main.cpp
    include/readpaf/paf.c
    include/readpaf/sdict.c
    src/polisher.cpp
    src/overlap.cpp
    src/sequence.cpp
    src/Sequence.cpp
    src/window.cpp
    src/main/preqclr.cpp)

target_include_directories(preqclrhtzy PRIVATE include/zstr/src)
target_include_directories(preqclrhtzy PRIVATE include/readpaf)
target_include_directories(preqclrhtzy PRIVATE include/)
target_include_directories(preqclrhtzy PRIVATE src/)
target_include_directories(preqclrhtzy PRIVATE include/rapidjson/include/)

#SUBDIRS := include src include/readpaf include/rapidjson/include include/zstr/src include/spoa

#add_subdirectory(vendor/rapidjson EXCLUDE_FROM_ALL)
add_subdirectory(vendor/bioparser EXCLUDE_FROM_ALL)
add_subdirectory(vendor/spoa EXCLUDE_FROM_ALL)
add_subdirectory(vendor/thread_pool EXCLUDE_FROM_ALL)
add_subdirectory(vendor/edlib EXCLUDE_FROM_ALL)

target_link_libraries(preqclrhtzy bioparser spoa thread_pool pthread edlib_static)

install(TARGETS preqclrhtzy DESTINATION bin)

if (racon_build_tests)
    set(racon_test_data_path ${PROJECT_SOURCE_DIR}/test/data/)
    configure_file("${PROJECT_SOURCE_DIR}/test/racon_test_config.h.in"
        "${PROJECT_BINARY_DIR}/config/racon_test_config.h")
    include_directories(${PROJECT_BINARY_DIR}/config)
    include_directories(${PROJECT_SOURCE_DIR}/src)

    add_executable(racon_test
        test/racon_test.cpp
        src/polisher.cpp
        src/overlap.cpp
        src/sequence.cpp
        src/window.cpp)

    add_subdirectory(vendor/googletest/googletest EXCLUDE_FROM_ALL)

    target_link_libraries(racon_test bioparser spoa thread_pool pthread
        edlib_static gtest_main)
endif(racon_build_tests)

if (racon_build_wrapper)
    set(racon_path ${PROJECT_BINARY_DIR}/bin/racon)
    set(rampler_path ${PROJECT_BINARY_DIR}/vendor/rampler/bin/rampler)
    configure_file(${PROJECT_SOURCE_DIR}/scripts/racon_wrapper.py
        ${PROJECT_BINARY_DIR}/${CMAKE_FILES_DIRECTORY}/racon_wrapper)
    file(COPY ${PROJECT_BINARY_DIR}/${CMAKE_FILES_DIRECTORY}/racon_wrapper
        DESTINATION ${PROJECT_BINARY_DIR}/bin
        FILE_PERMISSIONS OWNER_READ OWNER_EXECUTE GROUP_READ GROUP_EXECUTE
        WORLD_READ WORLD_EXECUTE)

    add_subdirectory(vendor/rampler)
endif(racon_build_wrapper)
