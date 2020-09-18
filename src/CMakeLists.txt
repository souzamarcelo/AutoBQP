cmake_minimum_required(VERSION 2.8)
include (InstallRequiredSystemLibraries)

project(autobqp)

set(CMAKE_BUILD_TYPE Release)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wno-unknown-pragmas -fPIC -DIL_STD -fexceptions -std=c++0x -O3 -DNDEBUG")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -g -ggdb")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -march=native -mtune=native -openmp -O3 -DNDEBUG")

find_package(Boost 1.50.0 COMPONENTS program_options timer system chrono filesystem)
include_directories(SYSTEM ${Boost_INCLUDE_DIRS})

##set(Boost_NO_SYSTEM_PATHS OFF)
#set(Boost_NO_SYSTEM_PATHS ON)
#set(Boost_USE_STATIC_LIBS ON)
#set(BOOST_INCLUDEDIR /usr/local/include)
#include_directories(SYSTEM ${BOOST_INCLUDEDIR})
#find_package(Boost COMPONENTS program_options timer system chrono filesystem)
#include_directories("${PROJECT_BINARY_DIR}")

add_executable(grammar grammar.cpp)
add_executable(runner runner.cpp)
target_link_libraries(grammar ${Boost_LIBRARIES} rt)
target_link_libraries(runner ${Boost_LIBRARIES} rt)