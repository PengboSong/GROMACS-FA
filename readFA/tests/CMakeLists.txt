cmake_minimum_required(VERSION 3.1.0)

add_executable(rvecTest rvecTest.cpp)

target_link_libraries(rvecTest FAcore)

target_link_libraries(rvecTest FAclib)

add_executable(containerTest containerTest.cpp)

target_link_libraries(containerTest FAcore)

add_executable(coordTest coordTest.cpp)

target_link_libraries(coordTest FAcore)

target_link_libraries(coordTest FAclib)
