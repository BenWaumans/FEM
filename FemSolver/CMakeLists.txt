list(REMOVE_DUPLICATES LIBRARIES)
list(REMOVE_DUPLICATES INCLUDE_DIRS)
message(STATUS "LIBRARIES = ${LIBRARIES}")
message(STATUS "INCLUDE_DIRS = ${INCLUDE_DIRS}")
include_directories(${INCLUDE_DIRS})

file(GLOB SOURCES "*.cpp")
message(STATUS "SOURCES: ${SOURCES}")

if(USE_CUDA)
  foreach(file IN LISTS SOURCES)
    set_property(SOURCE ${file} PROPERTY LANGUAGE CUDA)
  endforeach()
endif()

# configure a header file to pass some of the CMake settings
# to the source code
configure_file (
  "${PROJECT_SOURCE_DIR}/FemSolver/config/Config.hpp.in"
  "${PROJECT_BINARY_DIR}/FemSolver/Config.hpp"
  )

include_directories("${PROJECT_BINARY_DIR}/FemSolver")

# Add one executable per cpp file
add_executable(fem_solver ${SOURCES})
add_dependencies(fem_solver mfem_project opencascade_project)
target_link_libraries(fem_solver ${LIBRARIES})

install(TARGETS fem_solver DESTINATION ${CMAKE_INSTALL_PREFIX}/bin)