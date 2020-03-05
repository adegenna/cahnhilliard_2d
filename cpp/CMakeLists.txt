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

# Set paths to external library files here
set( MPI_LIBRARIES "/home/adegennaro/Projects/AEOLUS/cahnhilliard_2d/cpp/build/external/petsc/arch-linux2-c-opt/lib/libmpi.so" )
message( STATUS "MPI_LIBRARIES ${MPI_LIBRARIES}")
set( HDF5_LIBRARIES "/home/adegennaro/Projects/AEOLUS/cahnhilliard_2d/cpp/build/external/petsc/arch-linux2-c-opt/lib/libhdf5.so" )
message( STATUS "HDF5_LIBRARIES ${HDF5_LIBRARIES}")

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

if (NOT CMAKE_BUILD_TYPE)
  # set(CMAKE_BUILD_TYPE "Release")
  set(CMAKE_BUILD_TYPE "Debug")
endif()

add_definitions( -DSRCDIR=\"${CMAKE_SOURCE_DIR}/\")

set(CMAKE_CXX_FLAGS_DEBUG "-O0 -g -Wall -Werror=reorder")
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -g -DBOOST_LOG_DYN_LINK -DNDEBUG")

set( EXTERNAL_INSTALL_DIR ${CMAKE_BINARY_DIR}/external_install )
add_subdirectory(external)

message( STATUS "PETSC_INCLUDES ${PETSC_INCLUDES}")
message( STATUS "PETSC_LIBRARIES ${PETSC_LIBRARIES}")

include_directories( ${CMAKE_SOURCE_DIR}
                     ${Boost_INCLUDE_DIRS}
		     ${PETSC_INCLUDES}
	     
)

#############################################
# Build openmp solver
#############################################

add_library( ch_src
	     src/stochastic_euler.hpp
	     src/chparams.cpp
	     src/utils_ch.cpp
	     src/cahnhilliard.cpp
	     src/cahnhilliard_thermal.cpp
	     src/cahnhilliard_thermal_nodiffusion.cpp
	     src/run_ch_solver.cpp
	     src/right_hand_side.cpp  
  )

add_executable( ch2d
  src/driver.cpp
  )

target_link_libraries( ch2d
  ${Boost_LIBRARIES}   
  ch_src
)

#############################################
# Swig/Python build
#############################################

FIND_PACKAGE( SWIG REQUIRED )
INCLUDE( ${SWIG_USE_FILE} )

FIND_PACKAGE( PythonLibs )
INCLUDE_DIRECTORIES( ${PYTHON_INCLUDE_PATH} )
INCLUDE_DIRECTORIES( ${CMAKE_CURRENT_SOURCE_DIR} )
INCLUDE_DIRECTORIES( ${CMAKE_CURRENT_SOURCE_DIR}/src )
	    
SET( CMAKE_SWIG_FLAGS "-threads" )
SET( CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS} -fPIC")

SET_SOURCE_FILES_PROPERTIES( src/cahnhilliard.i PROPERTIES CPLUSPLUS ON )
SET_SOURCE_FILES_PROPERTIES( src/cahnhilliard.i PROPERTIES SWIG_FLAGS "-threads" )

#SWIG_ADD_LIBRARY( cahnhilliard TYPE SHARED LANGUAGE python SOURCES src/cahnhilliard.i )
SWIG_ADD_MODULE( cahnhilliard python src/cahnhilliard.i src/chparams.cpp src/utils_ch.cpp )

SWIG_LINK_LIBRARIES( cahnhilliard
                     ${PYTHON_LIBRARIES}
		     ch_src
)
		   
# # Copy the shared objects and Swig wrappers to cpp/swig/
# add_custom_target( mvcmd_so ALL
#   COMMAND ${CMAKE_COMMAND} -E copy "${CMAKE_SOURCE_DIR}/build/_cahnhilliard.so" "${CMAKE_SOURCE_DIR}/swig/_cahnhilliard.so"
#   DEPENDS "${CMAKE_SOURCE_DIR}/build/_cahnhilliard.so"
# )
# add_custom_target( mvcmd_py ALL
#   COMMAND ${CMAKE_COMMAND} -E copy "${CMAKE_SOURCE_DIR}/build/cahnhilliard.py" "${CMAKE_SOURCE_DIR}/swig/cahnhilliard.py"
#   DEPENDS "${CMAKE_SOURCE_DIR}/build/cahnhilliard.py"
# )


#############################################
# Petsc-3D split solver
#############################################

add_executable( ch3d_split
  implicit_petsc_3d/ch3d_split.cpp
  implicit_petsc_3d/utils_ch_implicit.cpp
  implicit_petsc_3d/boundary_conditions.cpp
  implicit_petsc_3d/initial_conditions.cpp
  implicit_petsc_3d/rhs_implicit.cpp
  implicit_petsc_3d/petsc_event_handling.cpp
  implicit_petsc_3d/temperature_dependence.cpp
)

add_dependencies(ch3d_split petsc_external)

target_link_libraries( ch3d_split
  m
  ${MPI_LIBRARIES}
  ${HDF5_LIBRARIES}
  ${PETSC_LIBRARIES}
)

#############################################
# Postprocessing
#############################################

add_executable( postprocess
  implicit_petsc/postprocess.cpp
  implicit_petsc/utils_ch_implicit.cpp
)

add_dependencies(postprocess petsc_external)

target_link_libraries( postprocess
  m
  ${MPI_LIBRARIES}
  ${HDF5_LIBRARIES}
  ${PETSC_LIBRARIES}
)

#############################################
# Preprocessing
#############################################

add_executable( preprocess
  implicit_petsc/preprocess.cpp
  implicit_petsc/utils_ch_implicit.cpp
)

add_dependencies( preprocess petsc_external )

target_link_libraries( preprocess
  m
  ${MPI_LIBRARIES}
  ${HDF5_LIBRARIES}
  ${PETSC_LIBRARIES}
)

#############################################
# Preprocessing-3d
#############################################

add_executable( preprocess_3d
  implicit_petsc_3d/preprocess.cpp
  implicit_petsc_3d/utils_ch_implicit.cpp
)

add_dependencies( preprocess_3d petsc_external )

target_link_libraries( preprocess_3d
  m
  ${MPI_LIBRARIES}
  ${HDF5_LIBRARIES}
  ${PETSC_LIBRARIES}
)

# # Copy petsc implicit solver driver executable to implicit_petsc/examples/
# add_custom_target( mvcmd_implicit ALL
#   COMMAND ${CMAKE_COMMAND} -E copy "${CMAKE_SOURCE_DIR}/build/jfnk_2d" "${CMAKE_SOURCE_DIR}/implicit_petsc/examples/jfnk_2d"
#   DEPENDS "${CMAKE_SOURCE_DIR}/build/jfnk_2d"
#   )

#############################################
# PETSc ex15
#############################################

add_executable( ex15
  implicit_petsc/ex15.cpp
)

add_dependencies( ex15 petsc_external )

target_link_libraries( ex15
  m
  ${MPI_LIBRARIES}
  ${HDF5_LIBRARIES}
  ${PETSC_LIBRARIES}
)



#############################################
# Petsc-2D implicit solver
#############################################

add_executable( ch2d_split
  implicit_petsc/ch2d_split.cpp
  implicit_petsc/utils_ch_implicit.cpp
  implicit_petsc/boundary_conditions.cpp
  implicit_petsc/initial_conditions.cpp
  implicit_petsc/rhs_implicit.cpp
  implicit_petsc/petsc_event_handling.cpp
  implicit_petsc/temperature_dependence.cpp
)

add_dependencies(ch2d_split petsc_external)

target_link_libraries( ch2d_split
  m
  ${MPI_LIBRARIES}
  ${HDF5_LIBRARIES}
  ${PETSC_LIBRARIES}
)


#############################################
# Petsc-2D explicit solver
#############################################

add_executable( ch2d_explicit
  explicit_petsc/ch2d_explicit.cpp
  explicit_petsc/utils_ch_explicit.cpp
  explicit_petsc/boundary_conditions.cpp
  explicit_petsc/initial_conditions.cpp
  explicit_petsc/rhs_explicit.cpp
  explicit_petsc/petsc_event_handling.cpp
  explicit_petsc/temperature_dependence.cpp
)

add_dependencies(ch2d_explicit petsc_external)

target_link_libraries( ch2d_explicit
  m
  ${MPI_LIBRARIES}
  ${HDF5_LIBRARIES}
  ${PETSC_LIBRARIES}
)

#############################################
# MMS test: Petsc-2D explicit solver 
#############################################

add_executable( mmstest_ch2d_explicit
  explicit_petsc/mmstest_ch2d_explicit.cpp
  explicit_petsc/utils_ch_explicit.cpp
  explicit_petsc/boundary_conditions.cpp
  explicit_petsc/initial_conditions.cpp
  explicit_petsc/rhs_explicit.cpp
  explicit_petsc/petsc_event_handling.cpp
  explicit_petsc/temperature_dependence.cpp
)

add_dependencies(mmstest_ch2d_explicit petsc_external)

target_link_libraries( mmstest_ch2d_explicit
  m
  ${MPI_LIBRARIES}
  ${HDF5_LIBRARIES}
  ${PETSC_LIBRARIES}
)
