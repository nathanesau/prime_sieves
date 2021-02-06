cmake_minimum_required(VERSION 3.10)

project(main)

set(CMAKE_CXX_STANDARD 17)
set(THREADS_PREFER_PTHREAD_FLAG ON)
set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/CMake/Modules)

find_package(GMP REQUIRED)
find_package(Threads REQUIRED)
find_package(Boost 1.75 COMPONENTS program_options REQUIRED)

INCLUDE_DIRECTORIES("${GMP_INCLUDES}")
INCLUDE_DIRECTORIES(${Boost_INCLUDE_DIR})

add_executable(main
    src/main.cpp
    src/params.h
    src/sieve.h
)

target_link_libraries(main ${GMP_LIBRARIES} ${Boost_LIBRARIES} Threads::Threads m)
