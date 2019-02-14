cmake_minimum_required(VERSION 3.6)
project(RELI)

set(CMAKE_CXX_STANDARD 11)
find_package(GSL REQUIRED)
include_directories (include)

set(SOURCE_FILES
    src/RELI.cpp
    src/RELI_impl.cpp
        include/RELI_impl.h)

add_executable(RELI ${SOURCE_FILES})
target_link_libraries(RELI GSL::gsl)