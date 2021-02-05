# https://github.com/mrsarm/gmp_test/blob/master/CMakeLists.txt
cmake_minimum_required(VERSION 3.1)
project(main C)

set(CMAKE_C_STANDARD 99)
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS}")
set(THREADS_PREFER_PTHREAD_FLAG ON)

find_package(Threads REQUIRED)

INCLUDE_DIRECTORIES(
    "${CMAKE_SOURCE_DIR}/src"
)

add_executable(main
    src/main.c
)

target_link_libraries(main PRIVATE Threads::Threads)
