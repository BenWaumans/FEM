cmake_minimum_required (VERSION 3.10)
PROJECT(fem_solver C CXX)

# set(MFEM_CMAKE_MODULE_PATH "mfem" CACHE PATH "Path to the mfem cmake module folder")


set(FEM_SOLVER_VERSION_MAJOR 0)
set(FEM_SOLVER_VERSION_MINOR 0)
set(FEM_SOLVER_VERSION_PATCH 0)

set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/config/cmake/modules)

find_package(MPI REQUIRED)
set(MPI_CXX_INCLUDE_DIRS ${MPI_CXX_INCLUDE_PATH})
find_package(HYPRE REQUIRED)
find_package(METIS REQUIRED)
find_package(OpenMP REQUIRED)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()
set(MFEM_DIR "mfem" CACHE PATH "Install dir of mfem")
# find_package(MFEM REQUIRED)
# find_package(mfem REQUIRED NAMES MFEM)
find_package(MFEM REQUIRED NAMES MFEM HINTS "${MFEM_DIR}"
                "${MFEM_DIR}/lib/cmake/mfem" NO_DEFAULT_PATH)

# List all possible libraries in order of dependencies.
# [METIS < SuiteSparse]:
#    With newer versions of SuiteSparse which include METIS header using 64-bit
#    integers, the METIS header (with 32-bit indices, as used by mfem) needs to
#    be before SuiteSparse.
set(FEM_SOLVER_TPLS MPI_CXX OpenMP METIS HYPRE MFEM)
# Add all *_FOUND libraries in the variable TPL_LIBRARIES.
set(TPL_LIBRARIES "")
set(TPL_INCLUDE_DIRS "")
foreach(TPL IN LISTS FEM_SOLVER_TPLS)
  if (${TPL}_FOUND)
    message(STATUS "FEM_SOLVER: using package ${TPL}")
    list(APPEND TPL_LIBRARIES ${${TPL}_LIBRARIES})
    list(APPEND TPL_INCLUDE_DIRS ${${TPL}_INCLUDE_DIRS})
  endif()
endforeach(TPL)
list(REMOVE_DUPLICATES TPL_LIBRARIES)
list(REMOVE_DUPLICATES TPL_INCLUDE_DIRS)
message(STATUS "TPL_LIBRARIES = ${TPL_LIBRARIES}")
message(STATUS "TPL_INCLUDE_DIRS = ${TPL_INCLUDE_DIRS}")
include_directories(${TPL_INCLUDE_DIRS})

get_property(dirs DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} PROPERTY INCLUDE_DIRECTORIES)
foreach(dir ${dirs})
  message(STATUS "dir='${dir}'")
endforeach()

# configure a header file to pass some of the CMake settings
# to the source code
configure_file (
  "${PROJECT_SOURCE_DIR}/Config.h.in"
  "${PROJECT_BINARY_DIR}/Config.h"
  )

# Add one executable per cpp file
add_executable(fem_solver test1.cpp)
target_link_libraries(fem_solver PUBLIC ${TPL_LIBRARIES})


# add the binary tree to the search path for include files
# so that we will find TutorialConfig.h
include_directories("${PROJECT_BINARY_DIR}")
install(TARGETS fem_solver DESTINATION ${CMAKE_INSTALL_PREFIX}/bin)