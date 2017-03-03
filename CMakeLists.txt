set(CMAKE_VERBOSE_MAKEFILE 1)
cmake_minimum_required(VERSION 2.6)

# Set a default build type if none was specified
if(NOT CMAKE_BUILD_TYPE)
  message(STATUS "Setting build type to 'Release' as none was specified.")
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build." FORCE)
endif()
message("Building in ${CMAKE_BUILD_TYPE} mode")

set(CMAKE_CXX_FLAGS "--std=c++14")
set(CMAKE_CXX_FLAGS_DEBUG "-O0 -ggdb -g -p")
set(CMAKE_CXX_FLAGS_RELEASE "-ggdb -Ofast -fstrict-aliasing -DNDEBUG -march=native")
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-g -ggdb -Ofast -fstrict-aliasing -march=native")

project(OnlineRLBWT)

add_subdirectory(${PROJECT_SOURCE_DIR}/extern/Basics)
include_directories(${PROJECT_SOURCE_DIR}/extern/Basics)
add_executable(${PROJECT_NAME} OnlineRLBWT.cpp DynRLE.cpp)
target_link_libraries(${PROJECT_NAME} Basics)
