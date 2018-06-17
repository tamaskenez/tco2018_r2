cmake_minimum_required(VERSION 3.2)
project(tco2018_2)

set(CMAKE_CXX_STANDARD 14)
add_executable(fortester fortester.cpp CrystalLighting.h)
add_executable(testlast testlast.cpp CrystalLighting.h)
