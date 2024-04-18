cmake_minimum_required(VERSION 3.1.0)

file(GLOB READFA_SRC *.cpp)

add_library(FAcore STATIC ${READFA_SRC})
