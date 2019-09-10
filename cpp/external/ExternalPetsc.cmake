#########################################################
message( STATUS "Building external Petsc project.")
#########################################################

if(${CMAKE_BUILD_TYPE} STREQUAL "Release")
    set(PETSC_OPT_FLAGS "-O3")
    set(PETSC_DEBUGGING "no")
else()
    set(PETSC_OPT_FLAGS "-O0")
    set(PETSC_DEBUGGING "yes")
endif()

find_package( PythonInterp 2.7 REQUIRED )

ExternalProject_Add(
        petsc_external

	URL http://ftp.mcs.anl.gov/pub/petsc/release-snapshots/petsc-lite-3.9.4.tar.gz

	#URL http://ftp.mcs.anl.gov/pub/petsc/release-snapshots/petsc-lite-3.8.3.tar.gz
        #URL_MD5 b21d15cdc033d50cc47b9997ef490fef
        
        BUILD_IN_SOURCE 1
        SOURCE_DIR ${CMAKE_BINARY_DIR}/external/petsc/

        CONFIGURE_COMMAND
        ${PYTHON_EXECUTABLE} ${CMAKE_BINARY_DIR}/external/petsc/configure
        PETSC_DIR=${CMAKE_BINARY_DIR}/external/petsc
        PETSC_ARCH=arch-linux2-c-opt
	--with-fc=0 --with-pic=1 --download-metis MAKEFLAGS=$MAKEFLAGS  COPTFLAGS=${PETSC_OPT_FLAGS} CXXOPTFLAGS=${PETSC_OPT_FLAGS} --with-debugging=${PETSC_DEBUGGING} ${PETSC_64_BIT_INDEX_FLAG} --download-openmpi --download-hdf5
        #--with-cc=${MPI_C_COMPILER} --with-cxx=${MPI_CXX_COMPILER} --with-fc=0 --with-pic=1 --download-metis MAKEFLAGS=$MAKEFLAGS  COPTFLAGS=${PETSC_OPT_FLAGS} CXXOPTFLAGS=${PETSC_OPT_FLAGS} --with-debugging=${PETSC_DEBUGGING} ${PETSC_64_BIT_INDEX_FLAG} --download-openmpi

        BUILD_COMMAND
        make
        PETSC_DIR=${CMAKE_BINARY_DIR}/external/petsc/
        PETSC_ARCH=arch-linux2-c-opt all

        INSTALL_COMMAND ""
)

ExternalProject_Get_Property(petsc_external source_dir)

set(PETSC_INCLUDES
        ${source_dir}/include
        ${source_dir}/arch-linux2-c-opt/include
        CACHE PATH "")

set(PETSC_LIBRARIES
        ${source_dir}/arch-linux2-c-opt/lib/libpetsc.so
        CACHE FILEPATH "")

set(METIS_INCLUDE_DIR ${source_dir}/arch-linux2-c-opt/include CACHE PATH "")
set(METIS_LIBRARY ${source_dir}/arch-linux2-c-opt/lib/libpetsc.so CACHE FILEPATH "")
