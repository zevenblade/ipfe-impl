cmake_minimum_required(VERSION 3.11)
project(cifer C)

# C compiler settings
set(CMAKE_C_STANDARD 11)

set(CMAKE_C_FLAGS "-Wall -Wextra -Wpedantic -Wredundant-decls \
-Wnested-externs -Wmissing-include-dirs -Wstrict-prototypes -Wwrite-strings \
-Wno-unused-parameter -Wshadow -Wold-style-definition")

# GCC warnings that Clang doesn't provide
if ("${CMAKE_C_COMPILER_ID}" STREQUAL "GNU")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wjump-misses-init -Wlogical-op")
endif ()

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif ()

set(CMAKE_C_FLAGS_DEBUG "-g")
set(CMAKE_C_FLAGS_RELEASE "-O2")

# Header files
include_directories(
        include
        external/munit
        external/uthash
        external/
        )

link_directories(external/amcl/)

# Library sources
set(library_SOURCES
        src/data/mat.c
        src/data/mat_curve.c
        src/data/vec.c
        src/data/vec_float.c
        src/data/vec_curve.c
        src/internal/big.c
        src/internal/common.c
        src/internal/dlog.c
        src/internal/hash.c
        src/internal/keygen.c
        src/internal/prime.c
        src/internal/str.c
        src/innerprod/fullysec/uzpipfe.c
        src/sample/normal.c
        src/sample/normal_cumulative.c
        src/sample/normal_double.c
        src/sample/normal_double_constant.c
        src/sample/normal_cdt.c
        src/sample/normal_negative.c
        src/sample/uniform.c
        src/abe/policy.c
        src/abe/uabipfe.c
        )

add_library(cifer SHARED ${library_SOURCES})

# Link libraries that are used in our library
target_link_libraries(cifer gmp sodium m amcl)
# Search for protobuf-c library
find_package(PkgConfig)
if (PKG_CONFIG_FOUND)
    pkg_check_modules(PBUFC "libprotobuf-c >= 1.0.0")
    if (PBUFC_FOUND)
        target_link_libraries(cifer ${PBUFC_LIBRARIES})
        add_definitions(${PBUFC_CFLAGS} ${PBUFC_CFLAGS_OTHER})
        link_directories(${PBUFC_LIBRARY_DIRS})
    endif()
endif()

# Install library and copy header to install dir
install(TARGETS cifer DESTINATION lib)
install(DIRECTORY include/ DESTINATION include)

# Create an executable
set(binary_SOURCES
        )

#add_executable(cifer_test ${binary_SOURCES})

#target_link_libraries(cifer_test PRIVATE cifer)

#add_custom_target(test COMMAND cifer_test VERBATIM)

add_custom_target(docs COMMAND doxygen WORKING_DIRECTORY .. VERBATIM)
