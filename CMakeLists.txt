cmake_minimum_required (VERSION 3.8 FATAL_ERROR)
PROJECT(FEM_SOLVER VERSION 0.0.0.0 LANGUAGES C CXX)

set(EP_BASE "")
set_directory_properties(PROPERTIES EP_BASE "${CMAKE_BINARY_DIR}/externals")

# Require C++11 and disable compiler-specific extensions
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

option(USE_MPI "Enable MPI" OFF)
if(USE_MPI)
  option(BUILD_MPI "Build MPI" OFF)
endif()

option(USE_CUDA "Enable CUDA" OFF)
if (USE_CUDA)
  set(CUDA_ARCH
    "sm_30" "sm_32" "sm_35" "sm_37"
    "sm_50" "sm_52" "sm_53"
    "sm_60" "sm_61" "sm_62"
    "sm_70" "sm_72" "sm_75"
    CACHE STRING "Target CUDA architecture(s).")
  enable_language(CUDA)
  message(STATUS "Using CUDA architecture(s): ${CUDA_ARCH}")

  set(CMAKE_CUDA_STANDARD 11)
  set(CMAKE_CUDA_STANDARD_REQUIRED ON)
  set(CMAKE_CUDA_EXTENSIONS OFF)
  foreach(arch ${CUDA_ARCH})
    set(APPEND CUDA_ARCHS "-arch=${arch}")
  endforeach()

  set(CMAKE_CUDA_FLAGS ${CUDA_ARCHS} "--expt-extended-lambda"
      CACHE STRING "CUDA flags set for FEM_SOLVER" FORCE)
  if(USE_MPI)
    set(CUDA_CCBIN_COMPILER ${MPI_CXX_COMPILER})
  else()
    set(CUDA_CCBIN_COMPILER ${CMAKE_CXX_COMPILER})
  endif()
  if (NOT MSVC)
    string(APPEND CMAKE_CUDA_FLAGS " -ccbin ${CUDA_CCBIN_COMPILER}")
  endif()
endif()

find_package(OpenMP REQUIRED)
if (OPENMP_FOUND)
  message(STATUS "FEM_SOLVER: using package OpenMP")
  set(LIBRARIES ${OPENMP_LIBRARIES})
  set(INCLUDE_DIRS ${OPENMP_INCLUDE_DIRS})
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
  set(USE_OPENMP ON)
else()
  set(USE_OPENMP OFF)
endif()

set(DEPENDENCIES_INSTALL_DIR "${PROJECT_BINARY_DIR}/dependencies")
set(DEPENDENCIES_LIB_DIR "${DEPENDENCIES_INSTALL_DIR}/lib")
set(DEPENDENCIES_BIN_DIR "${DEPENDENCIES_INSTALL_DIR}/bin")
set(DEPENDENCIES_INCLUDE_DIR "${DEPENDENCIES_INSTALL_DIR}/include")

macro(add_static_lib TARGET LIB)
  ExternalProject_Get_Property(${TARGET} INSTALL_DIR)
  if(NOT "${LIB}" IN_LIST LIBRARIES)
    add_library(${LIB} STATIC IMPORTED)
    set_target_properties(${LIB} PROPERTIES IMPORTED_LOCATION
      "${INSTALL_DIR}/lib/${CMAKE_STATIC_LIBRARY_PREFIX}${LIB}${CMAKE_STATIC_LIBRARY_SUFFIX}")
    list(APPEND LIBRARIES "${LIB}")
  endif()
endmacro(add_static_lib)

macro(add_include TARGET)
  ExternalProject_Get_Property(${TARGET} INSTALL_DIR)
  set(SUB_DIR "${ARGV1}")
  if(SUB_DIR)
    list(APPEND INCLUDE_DIRS "${INSTALL_DIR}/include/${SUB_DIR}")
  else()
    list(APPEND INCLUDE_DIRS "${INSTALL_DIR}/include")
  endif()
  set(${TARGET}_INSTALL_DIR "${INSTALL_DIR}")
endmacro(add_include)

# Set use compiler flags for configure based projects
set(C_FLAGS ${CMAKE_C_FLAGS} ${CMAKE_C_FLAGS_RELEASE})
set(CXX_FLAGS ${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_RELEASE})
set(EXE_LINKER_FLAGS ${CMAKE_EXE_LINKER} ${CMAKE_EXE_LINKER})

set(FEM_ARGS "-DCMAKE_INSTALL_PREFIX=<INSTALL_DIR>"
             "-DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}"
             "-DCMAKE_AR=${CMAKE_AR}"
             "-DCMAKE_LINKER=${CMAKE_LINKER}"
             "-DCMAKE_NM=${CMAKE_NM}"
             "-DCMAKE_OBJCOPY=${CMAKE_OBJCOPY}"
             "-DCMAKE_OBJDUMP=${CMAKE_OBJDUMP}"
             "-DCMAKE_RANLIB=${CMAKE_RANLIB}"
             "-DCMAKE_STRIP=${CMAKE_STRIP}"
             "-DCMAKE_EXE_LINKER_FLAGS=${CMAKE_EXE_LINKER_FLAGS}"
             "-DCMAKE_EXE_LINKER_FLAGS_DEBUG=${CMAKE_EXE_LINKER_FLAGS_DEBUG}"
             "-DCMAKE_EXE_LINKER_FLAGS_MINSIZEREL=${CMAKE_EXE_LINKER_FLAGS_MINSIZEREL}"
             "-DCMAKE_EXE_LINKER_FLAGS_RELEASE=${CMAKE_EXE_LINKER_FLAGS_RELEASE}"
             "-DCMAKE_EXE_LINKER_FLAGS_RELWITHDEBINFO=${CMAKE_EXE_LINKER_FLAGS_RELWITHDEBINFO}"
             "-DCMAKE_MODULE_LINKER_FLAGS=${CMAKE_MODULE_LINKER_FLAGS}"
             "-DCMAKE_MODULE_LINKER_FLAGS_DEBUG=${CMAKE_MODULE_LINKER_FLAGS_DEBUG}"
             "-DCMAKE_MODULE_LINKER_FLAGS_MINSIZEREL=${CMAKE_MODULE_LINKER_FLAGS_MINSIZEREL}"
             "-DCMAKE_MODULE_LINKER_FLAGS_RELEASE=${CMAKE_MODULE_LINKER_FLAGS_RELEASE}"
             "-DCMAKE_MODULE_LINKER_FLAGS_RELWITHDEBINFO=${CMAKE_MODULE_LINKER_FLAGS_RELWITHDEBINFO}"
             "-DCMAKE_SHARED_LINKER_FLAGS=${CMAKE_SHARED_LINKER_FLAGS}"
             "-DCMAKE_SHARED_LINKER_FLAGS_DEBUG=${CMAKE_SHARED_LINKER_FLAGS_DEBUG}"
             "-DCMAKE_SHARED_LINKER_FLAGS_MINSIZEREL=${CMAKE_SHARED_LINKER_FLAGS_MINSIZEREL}"
             "-DCMAKE_SHARED_LINKER_FLAGS_RELEASE=${CMAKE_SHARED_LINKER_FLAGS_RELEASE}"
             "-DCMAKE_SHARED_LINKER_FLAGS_RELWITHDEBINFO=${CMAKE_SHARED_LINKER_FLAGS_RELWITHDEBINFO}"
             "-DCMAKE_STATIC_LINKER_FLAGS=${CMAKE_STATIC_LINKER_FLAGS}"
             "-DCMAKE_STATIC_LINKER_FLAGS_DEBUG=${CMAKE_STATIC_LINKER_FLAGS_DEBUG}"
             "-DCMAKE_STATIC_LINKER_FLAGS_MINSIZEREL=${CMAKE_STATIC_LINKER_FLAGS_MINSIZEREL}"
             "-DCMAKE_STATIC_LINKER_FLAGS_RELEASE=${CMAKE_STATIC_LINKER_FLAGS_RELEASE}"
             "-DCMAKE_STATIC_LINKER_FLAGS_RELWITHDEBINFO=${CMAKE_STATIC_LINKER_FLAGS_RELWITHDEBINFO}"
             "-DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}"
             "-DCMAKE_C_COMPILER_AR=${CMAKE_C_COMPILER_AR}"
             "-DCMAKE_C_COMPILER_RANLIB=${CMAKE_C_COMPILER_RANLIB}"
             "-DCMAKE_C_FLAGS=${CMAKE_C_FLAGS}"
             "-DCMAKE_C_FLAGS_DEBUG=${CMAKE_C_FLAGS_DEBUG}"
             "-DCMAKE_C_FLAGS_MINSIZEREL=${CMAKE_C_FLAGS_MINSIZEREL}"
             "-DCMAKE_C_FLAGS_RELEASE=${CMAKE_C_FLAGS_RELEASE}"
             "-DCMAKE_C_FLAGS_RELWITHDEBINFO=${CMAKE_C_FLAGS_RELWITHDEBINFO}"
             "-DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}"
             "-DCMAKE_CXX_COMPILER_AR=${CMAKE_CXX_COMPILER_AR}"
             "-DCMAKE_CXX_COMPILER_RANLIB=${CMAKE_CXX_COMPILER_RANLIB}"
             "-DCMAKE_CXX_FLAGS=${CMAKE_CXX_FLAGS}"
             "-DCMAKE_CXX_FLAGS_DEBUG=${CMAKE_CXX_FLAGS_DEBUG}"
             "-DCMAKE_CXX_FLAGS_MINSIZEREL=${CMAKE_CXX_FLAGS_MINSIZEREL}"
             "-DCMAKE_CXX_FLAGS_RELEASE=${CMAKE_CXX_FLAGS_RELEASE}"
             "-DCMAKE_CXX_FLAGS_RELWITHDEBINFO=${CMAKE_CXX_FLAGS_RELWITHDEBINFO}"
             "-DCMAKE_CXX_STANDARD=${CMAKE_CXX_STANDARD}")

include(ExternalProject)

if(USE_MPI)
  find_package(MPI REQUIRED)
  if (MPI_CXX_FOUND)
    list(APPEND LIBRARIES ${MPI_CXX_LIBRARIES})
    list(APPEND INCLUDE_DIRS ${MPI_CXX_INCLUDE_PATH})
  endif()

  set(MPI_COMPILER "-DMPIEXEC_EXECUTABLE=${MPIEXEC_EXECUTABLE}"
                   "-DMPIEXEC_MAX_NUMPROCS=${MPIEXEC_MAX_NUMPROCS}"
                   "-DMPIEXEC_NUMPROC_FLAG=${MPIEXEC_NUMPROC_FLAG}"
                   "-DMPIEXEC_POSTFLAGS=${MPIEXEC_POSTFLAGS}"
                   "-DMPIEXEC_PREFLAGS=${MPIEXEC_PREFLAGS}"
                   "-DMPI_CXX_ADDITIONAL_INCLUDE_DIRS=${MPI_CXX_ADDITIONAL_INCLUDE_DIRS}"
                   "-DMPI_CXX_COMPILER=${MPI_CXX_COMPILER}"
                   "-DMPI_CXX_COMPILE_DEFINITIONS=${MPI_CXX_COMPILE_DEFINITIONS}"
                   "-DMPI_CXX_COMPILE_OPTIONS=${MPI_CXX_COMPILE_OPTIONS}"
                   "-DMPI_CXX_HEADER_DIR=${MPI_CXX_HEADER_DIR}"
                   "-DMPI_CXX_LIB_NAMES=${MPI_CXX_LIB_NAMES}"
                   "-DMPI_CXX_LINK_FLAGS=${MPI_CXX_LINK_FLAGS}"
                   "-DMPI_CXX_SKIP_MPICXX=${MPI_CXX_SKIP_MPICXX}"
                   "-DMPI_C_ADDITIONAL_INCLUDE_DIRS=${MPI_C_ADDITIONAL_INCLUDE_DIRS}"
                   "-DMPI_C_COMPILER=${MPI_C_COMPILER}"
                   "-DMPI_C_COMPILE_DEFINITIONS=${MPI_C_COMPILE_DEFINITIONS}"
                   "-DMPI_C_COMPILE_OPTIONS=${MPI_C_COMPILE_OPTIONS}"
                   "-DMPI_C_HEADER_DIR=${MPI_C_HEADER_DIR}"
                   "-DMPI_C_LIB_NAMES=${MPI_C_LIB_NAMES}"
                   "-DMPI_C_LINK_FLAGS=${MPI_C_LINK_FLAGS}"
                   "-DMPI_mpi_LIBRARY=${MPI_mpi_LIBRARY}"
                   "-DMPI_mpi_cxx_LIBRARY=${MPI_mpi_cxx_LIBRARY}")


  set(HYPRE_ARGS  ${MPI_COMPILER}
                  "-DHYPRE_INSTALL_PREFIX=${DEPENDENCIES_INSTALL_DIR}"
                  "-DHYPRE_USING_OPENMP=${USE_OPENMP}"
                  "-S<SOURCE_DIR>/src")

  ExternalProject_Add(hypre_project
    URL "https://computation.llnl.gov/projects/hypre-scalable-linear-solvers-multigrid-methods/download/hypre-2.10.0b.tar.gz"
    CMAKE_ARGS ${FEM_ARGS} ${HYPRE_ARGS}
    # CONFIGURE_COMMAND ${CMAKE_COMMAND} ${FEM_ARGS} ${HYPRE_ARGS} <SOURCE_DIR>/src
  )

  add_static_lib(hypre_project HYPRE)
  add_include(hypre_project)

  set(METIS_ARGS ${MPI_COMPILER}
                 "-DGKLIB_PATH=<SOURCE_DIR>/GKlib"
                 "-DOPENMP=${USE_OPENMP}")

  ExternalProject_Add(metis_project
    URL "http://glaros.dtc.umn.edu/gkhome/fetch/sw/metis/metis-5.1.0.tar.gz"
    CMAKE_ARGS ${FEM_ARGS} ${METIS_ARGS}
  )

  add_static_lib(metis_project metis)
  add_include(metis_project)
endif()

set(MFEM_ARGS "-DMFEM_USE_CUDA=${USE_CUDA}"
              "-DMFEM_USE_MPI=${USE_MPI}"
              "-DMFEM_ENABLE_TESTING=ON"
              "-DMFEM_ENABLE_EXAMPLES=OFF"
              "-DMFEM_ENABLE_MINIAPPS=OFF"
              "-DMFEM_USE_OPENMP=${USE_OPENMP}"
              "-DMFEM_THREAD_SAFE=ON")
if (USE_MPI)
  list(APPEND MFEM_ARGS ${MPI_COMPILER}
                        "-DMETIS_DIR=${DEPENDENCIES_INSTALL_DIR}"
                        "-DHYPRE_DIR=${DEPENDENCIES_INSTALL_DIR}")
endif()
if (USE_CUDA)
  list(APPEND MFEM_ARGS "-DCUDA_ARCH=${CUDA_ARCH}")
endif()

if (USE_MPI)
  set(MFEM_DEPENDS hypre_project metis_project)
else()
  set(MFEM_DEPENDS "")
endif()

ExternalProject_Add(mfem_project
  # URL "https://github.com/mfem/mfem/archive/v4.0-rc2.tar.gz"
  GIT_REPOSITORY git@github.com:mfem/mfem.git
  GIT_TAG master
  GIT_SHALLOW ON
  GIT_PROGRESS ON
  CMAKE_ARGS ${FEM_ARGS} ${MFEM_ARGS}
  DEPENDS ${MFEM_DEPENDS}
)

add_static_lib(mfem_project mfem)
add_include(mfem_project)

list(APPEND OPENCASCADE_ARGS )
set(OPENCASCADE_ARGS "-DOCE_INSTALL_PREFIX=<INSTALL_DIR>"
                     "-DOCE_TESTING=OFF"
                     "-DOCE_VISUALISATION=OFF"
                     "-DOCE_ADD_HEADERS=ON"
                     "-DOCE_BUILD_SHARED_LIB=OFF"
                     "-DOCE_BUILD_STATIC_LIB=ON"
                     "-DOCE_DISABLE_X11=OFF")
if(USE_OPENMP)
  list(APPEND OPENCASCADE_ARGS "-DOCE_MULTITHREAD_LIBRARY=OPENMP")
endif()

ExternalProject_Add(opencascade_project
  URL "https://github.com/tpaviot/oce/archive/OCE-0.18.3.tar.gz"
  CMAKE_ARGS ${FEM_ARGS} ${OPENCASCADE_ARGS}
)

set(OPENCASCADE_LIBRARIES "TKernel"
                          "TKMath"
                          "TKG2d"
                          "TKG3d"
                          "TKGeomBase"
                          "TKBRep"
                          "TKGeomAlgo"
                          "TKTopAlgo"
                          "TKPrim"
                          "TKBO"
                          "TKHLR"
                          "TKMesh"
                          "TKShHealing"
                          "TKXMesh"
                          "TKBool"
                          "TKFillet"
                          "TKFeat"
                          "TKOffset"
                          "TKService"
                          "TKV3d"
                          "TKCDF"
                          "PTKernel"
                          "TKLCAF"
                          "FWOSPlugin"
                          "TKPShape"
                          "TKBinL"
                          "TKXmlL"
                          "TKPLCAF"
                          "TKTObj"
                          "TKShapeSchema"
                          "TKStdLSchema"
                          "TKCAF"
                          "TKBin"
                          "TKXml"
                          "TKPCAF"
                          "TKBinTObj"
                          "TKXmlTObj"
                          "TKStdSchema"
                          "TKSTL"
                          "TKXSBase"
                          "TKSTEPBase"
                          "TKIGES"
                          "TKSTEPAttr"
                          "TKSTEP209"
                          "TKSTEP"
                          "TKVRML"
                          "TKXCAF"
                          "TKXCAFSchema"
                          "TKXmlXCAF"
                          "TKBinXCAF"
                          "TKXDEIGES"
                          "TKXDESTEP")
list(REVERSE OPENCASCADE_LIBRARIES)
foreach(lib ${OPENCASCADE_LIBRARIES})
  add_static_lib(opencascade_project ${lib})
endforeach()
add_include(opencascade_project "oce")


set(NETGEN_ARGS "-DUSE_MPI=${USE_MPI}"
                "-DUSE_PYTHON=OFF"
                "-DUSE_GUI=OFF"
                "-DUSE_OCC=ON"
                "-DOCC_LIBRARY=${opencascade_project_INSTALL_DIR}/lib"
                "-DOCC_INCLUDE_DIR=${opencascade_project_INSTALL_DIR}/include/oce")

if (USE_MPI)
  set(NETGEN_DEPENDS opencascade_project metis_project)
else()
  set(NETGEN_DEPENDS opencascade_project)
endif()

ExternalProject_Add(netgen_project
  URL "https://github.com/NGSolve/netgen/archive/v6.2.1810.tar.gz"
  CMAKE_ARGS ${FEM_ARGS} ${NETGEN_ARGS}
  PATCH_COMMAND git -C <SOURCE_DIR> apply ${CMAKE_SOURCE_DIR}/patches/netgen.patch
  DEPENDS ${NETGEN_DEPENDS}
)

# add_static_lib(netgen_project netgen)
# add_include(netgen_project)

if (NOT USE_CUDA)
  list(APPEND GLVIS_ARGS "-DMFEM_DIR=${DEPENDENCIES_INSTALL_DIR}")

  ExternalProject_Add(glvis
    URL "http://glvis.github.io/releases/glvis-3.4.tgz"
    CMAKE_ARGS ${FEM_ARGS} ${GLVIS_ARGS}
    DEPENDS mfem_project
  )

  install(PROGRAMS "${DEPENDENCIES_BIN_DIR}/glvis${EXE_SUFFIX}" DESTINATION "${CMAKE_INSTALL_PREFIX}/bin")
endif()

add_subdirectory(FemSolver)
