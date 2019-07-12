cmake_minimum_required(VERSION 2.6)
project(ch2d)

include(ExternalProject)

find_package(Boost COMPONENTS program_options REQUIRED)
find_package(Threads)

find_package(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release")
endif()

add_definitions( -DSRCDIR=\"${CMAKE_SOURCE_DIR}/\")

set(CMAKE_CXX_FLAGS_DEBUG "-O0 -g -Wall -Werror=reorder")
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -g -DBOOST_LOG_DYN_LINK -DNDEBUG")


#############################################
# Build external dependencies
#############################################

set( EXTERNAL_INSTALL_DIR ${CMAKE_BINARY_DIR}/external_install )
add_subdirectory(external)

message( STATUS "PETSC_INCLUDES ${PETSC_INCLUDES}")

include_directories( ${CMAKE_SOURCE_DIR}
                     ${Boost_INCLUDE_DIRS}
                     ${PETSC_INCLUDES}
                     ${CMAKE_BINARY_DIR}/external/petsc/include/petsc/mpiuni # needed for non-mpi petsc version
)

add_executable( ch2d
  src/driver.cpp
  src/stochastic_euler.hpp
  src/chparams.cpp
  src/utils_ch.cpp
  src/cahnhilliard.cpp
  src/cahnhilliard_thermal.cpp
  src/cahnhilliard_thermal_nodiffusion.cpp
  src/run_ch_solver.cpp
  src/right_hand_side.cpp
)

add_dependencies(ch2d petsc_external)


target_link_libraries( ch2d
  ${Boost_LIBRARIES}   
  ${PETSC_LIBRARIES}
)
