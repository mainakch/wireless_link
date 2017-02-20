cmake_minimum_required(VERSION 3.7.2)
project (wireless_link)

set(PROJECT_SOURCE_DIR ${CMAKE_SOURCE_DIR}/src)
set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG}")

include_directories(include)
file(GLOB SOURCES "src/*.c")
find_library(GSL_LIB gsl)
find_library(BLAS_LIB cblas)
find_library(MATH_LIB m)

add_executable(wireless_link ${SOURCES})
target_link_libraries(wireless_link ${GSL_LIB} ${BLAS_LIB} ${MATH_LIB})
