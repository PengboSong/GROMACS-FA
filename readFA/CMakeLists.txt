cmake_minimum_required(VERSION 3.1.0)

project(readFA)

set(CMAKE_CXX_STANDARD   11)
set(CMAKE_CXX_STANDARD_REQUIRED   ON)

include_directories(${CMAKE_SOURCE_DIR}/include)

add_subdirectory(clib)

add_subdirectory(src)

add_subdirectory(tests)

add_subdirectory(tools)

add_executable(readFA readFA.cpp)

target_link_libraries(readFA FAcore)

target_link_libraries(readFA FAclib)

install(TARGETS readFA RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})

install(TARGETS readFMP RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})

install(TARGETS groupData RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
