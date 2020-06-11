cmake_minimum_required(VERSION 3.15.0)
project(TLC VERSION 0.1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 11)

add_executable(findInjective TLC.cpp)

set(NLopt_INCLUDE_DIR /usr/local/include CACHE PATH "NLopt include directory")
set(NLopt_LIB_DIR /usr/local/lib CACHE PATH "NLopt library directory")

target_include_directories(findInjective PUBLIC "${NLopt_INCLUDE_DIR}")

target_link_libraries(findInjective PUBLIC -L"${NLopt_LIB_DIR}" nlopt)

find_library(MATH_LIBRARY m)
if(MATH_LIBRARY)
    target_link_libraries(findInjective PUBLIC ${MATH_LIBRARY})
endif()

