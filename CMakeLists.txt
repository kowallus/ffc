cmake_minimum_required(VERSION 3.10)
project(FFC)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)

set(CMAKE_CXX_FLAGS_RELEASE  "${CMAKE_CXX_FLAGS} -O3 -Ofast -flto -funroll-loops -ffast-math")

if (NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

find_package(Threads REQUIRED)
link_libraries(Threads::Threads)

include_directories(CLI11/include)
include_directories(zstd)

set(HELPER_FILES
        utils/helper.cpp
        utils/helper.h)

set(SOURCES ${HELPER_FILES}
        ffc.cpp
        FASTA_Compress.cpp
        FASTA_Compress.h
        FASTA_Decompress.cpp
        FASTA_Decompress.h)

link_libraries(${CMAKE_CURRENT_SOURCE_DIR}/zstd/lib/libzstd.a)

add_executable(ffc ${SOURCES})
if(NOT CMAKE_SYSTEM_PROCESSOR MATCHES "^(arm|aarch64)")
    add_executable(ffc-noavx ${SOURCES})
    target_compile_definitions(ffc-noavx PUBLIC -DNO_AVX)
    target_compile_options(ffc PUBLIC -mavx2)
else()
    target_compile_definitions(ffc PUBLIC -DNO_AVX)
endif()
