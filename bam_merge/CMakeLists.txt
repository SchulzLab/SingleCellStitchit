cmake_minimum_required(VERSION 3.8)

project(bam_merge C CXX)

include_directories(include /usr/local/include/htslib/)
include_directories(include /usr/include/)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_INCLUDE_CURRENT_DIR ON)


add_compile_options(-Wall -pedantic -O3 -lstdc++fs)
add_compile_definitions(__ZLIB_AVAILABLE__)
add_compile_definitions(_FILE_OFFSET_BITS=64)
add_compile_definitions(__STDC_LIMIT_MACROS)

configure_file(CMakeLists.txt.in googletest-download/CMakeLists.txt)
execute_process(COMMAND "${CMAKE_COMMAND}" -G "${CMAKE_GENERATOR}" .
    WORKING_DIRECTORY "${CMAKE_BINARY_DIR}/googletest-download" )
execute_process(COMMAND "${CMAKE_COMMAND}" --build .
    WORKING_DIRECTORY "${CMAKE_BINARY_DIR}/googletest-download" )

add_subdirectory("${CMAKE_BINARY_DIR}/googletest-src"
                 "${CMAKE_BINARY_DIR}/googletest-build")

add_subdirectory(zlib)
add_subdirectory(bam)
add_subdirectory(core)
