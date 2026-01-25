cmake_minimum_required(VERSION 3.25)
project(chfilter)

set(CMAKE_CXX_STANDARD 17)

find_package(savvy REQUIRED)
find_package(Boost REQUIRED)

add_executable(chfilter main.cpp tcdf.hpp)
target_include_directories(chfilter PRIVATE include ${Boost_INCLUDE_DIRS})
target_link_libraries(chfilter savvy)