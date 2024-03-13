cmake_minimum_required(VERSION 3.23)
project(LamberBattin_C)

set(CMAKE_CXX_STANDARD 14)

include_directories(LamberBattin)

add_executable(LamberBattin_C
        LamberBattin/LAMBERBATTIN.h
        LamberBattin/LAMBERTBATTIN.cpp
        LamberBattin/main.cpp
        LamberBattin/seebatt.cpp
        LamberBattin/seebatt.h
        LamberBattin/seebattk.cpp
        LamberBattin/seebattk.h
        LamberBattin/vector.cpp
        LamberBattin/vector.h)
