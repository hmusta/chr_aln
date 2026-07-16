cmake_minimum_required(VERSION 3.16)
project(chr_aln)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release")
endif()
string(TOUPPER "${CMAKE_BUILD_TYPE}" CMAKE_BUILD_TYPE)

if(${CMAKE_BUILD_TYPE} MATCHES RELEASE)
    message(STATUS "Release mode")
elseif(${CMAKE_BUILD_TYPE} MATCHES RELWITHDEBINFO)
    message(STATUS "Release mode with debugging info")
    add_compile_options(-O3)
elseif(${CMAKE_BUILD_TYPE} MATCHES DEBUG)
    message(STATUS "Debug mode")
    add_compile_options(-Og -UNDEBUG -fno-omit-frame-pointer -g)
    add_link_options(-g)
elseif(${CMAKE_BUILD_TYPE} MATCHES PROFILE)
    set(PROFILING ON CACHE BOOL "Profile mode" FORCE)
    add_compile_options(-Og -DNDEBUG -fno-omit-frame-pointer -g -DPROFILING)
    add_link_options(-g)
endif()

add_compile_options(-Wall -march=native)

set(OPENMP ON CACHE BOOL "Enable OpenMP" FORCE)
find_package(OpenMP REQUIRED)

message(STATUS "Setting up UniAligner")
add_subdirectory(unialigner/tandem_aligner/src/tools)
add_library(tandem_aligner_lib STATIC
    unialigner/tandem_aligner/src/projects/tandem_aligner/min_interval.cpp
    unialigner/tandem_aligner/src/projects/tandem_aligner/lcp_interval.cpp
    unialigner/tandem_aligner/src/projects/tandem_aligner/cigar.cpp
)
set(TANDEM_DEP common sequence suffix_array_library min_queue_library)
set(TANDEM ${TANDEM_DEP} tandem_aligner_lib suffix_array min_queue)
foreach (X IN LISTS TANDEM)
    target_include_directories(${X} SYSTEM PRIVATE unialigner/tandem_aligner/src/tools)
    target_compile_features(${X} PRIVATE cxx_std_20)
    target_compile_options(${X} PRIVATE -lstdc++fs)
    target_link_options(${X} PRIVATE -static-libgcc -Wc++-compat -DHAVE_KALLOC -DKSW_CPU_DISPATCH -D_FILE_OFFSET_BITS=64 -fsigned-char)
endforeach()

target_link_libraries(tandem_aligner_lib PRIVATE ${TANDEM_DEP})

message(STATUS "Setting up WFA2-lib")
add_subdirectory(WFA2-lib)

message(STATUS "Setting up MUMmer")
add_library(mummer STATIC
    mummer/src/essaMEM/sparseSA.cpp
    mummer/src/essaMEM/sssort_compact.cc
)
target_include_directories(mummer SYSTEM PUBLIC mummer/include)

message(STATUS "Setting up C++ Progress Bar")
add_library(progress_bar STATIC
    cpp_progress_bar/progress_bar.cpp
)
target_include_directories(progress_bar SYSTEM PUBLIC cpp_progress_bar)
target_compile_features(progress_bar PRIVATE cxx_std_20)

message(STATUS "Setting up chr_aln")
add_library(libchraln
    helpers.cpp
    wfa_switch.cpp
    repeat_aligner.cpp
    chaining.cpp
)
target_compile_options(libchraln PRIVATE -Werror)
target_compile_features(libchraln PRIVATE cxx_std_20)
target_include_directories(libchraln PUBLIC ${CMAKE_SOURCE_DIR})
target_include_directories(libchraln SYSTEM PUBLIC
    WFA2-lib
    unialigner/tandem_aligner/src/tools
    unialigner/tandem_aligner/src/projects/tandem_aligner
    mummer/include/mummer
)
target_link_libraries(libchraln PRIVATE
    wfa2::wfa2cpp_static
    tandem_aligner_lib
    OpenMP::OpenMP_CXX
    progress_bar
    mummer
)

if(PROFILING)
    target_include_directories(libchraln SYSTEM PUBLIC
        gperftools/include
    )
    target_link_libraries(libchraln PRIVATE
        gperftools/lib/libprofiler.a
        libunwind/lib/libunwind.so
    )
endif()

add_executable(chr_aln
    main.cpp
)
target_compile_options(chr_aln PRIVATE -Werror)
target_compile_features(chr_aln PRIVATE cxx_std_20)
target_link_libraries(chr_aln PRIVATE
    libchraln
    OpenMP::OpenMP_CXX
    progress_bar
)

include(FetchContent)
FetchContent_Declare(
  googletest
  DOWNLOAD_EXTRACT_TIMESTAMP
  URL https://github.com/google/googletest/archive/refs/tags/v1.17.0.zip
)
# For Windows: Prevent overriding the parent project's compiler/linker settings
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(googletest)

enable_testing()

add_executable(unit_tests
    tests/test_wfa.cpp
    tests/test_chaining.cpp
    tests/test_seeding.cpp
)

target_compile_options(unit_tests PRIVATE -Werror)
target_compile_features(unit_tests PRIVATE cxx_std_20)
target_include_directories(unit_tests PRIVATE ${CMAKE_SOURCE_DIR})

target_link_libraries(unit_tests
    GTest::gtest_main
    libchraln
)

include(GoogleTest)
gtest_discover_tests(unit_tests)
