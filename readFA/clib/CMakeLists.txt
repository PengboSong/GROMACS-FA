cmake_minimum_required(VERSION 3.1.0)

file(GLOB READFA_CLIB_SRC *.c)

add_library(FAclib STATIC ${READFA_CLIB_SRC})
