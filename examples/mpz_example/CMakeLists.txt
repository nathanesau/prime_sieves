# https://github.com/mrsarm/gmp_test/blob/master/CMakeLists.txt
cmake_minimum_required(VERSION 3.1)
project(main C)

set(CMAKE_C_STANDARD 99)
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS}")

set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/CMake/Modules)
find_package(GMP REQUIRED)

INCLUDE_DIRECTORIES(
    "${GMP_INCLUDES}"
    "${CMAKE_SOURCE_DIR}/src"
)

add_executable(main
    src/main.c
)

target_link_libraries(main ${GMP_LIBRARIES})
