cmake_minimum_required(VERSION 3.0.2)

macro(use_cxx11)
  if (CMAKE_VERSION VERSION_LESS "3.1")
    if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
      set (CMAKE_CXX_FLAGS "--std=gnu++11 ${CMAKE_CXX_FLAGS}")
    endif ()
  else ()
    set (CMAKE_CXX_STANDARD 11)
  endif ()
endmacro(use_cxx11)

project (ACFcalculator)
use_cxx11()
add_executable(ACFcalculator ACFcalculator.cpp)
find_package(Boost 1.55 REQUIRED COMPONENTS program_options regex)
include_directories(${Boost_INCLUDE_DIRS})
target_link_libraries(ACFcalculator ${Boost_LIBRARIES})

OPTION (USE_OpenMP "Use OpenMP" ON)
IF(USE_OpenMP)
  FIND_PACKAGE(OpenMP)
  IF(OPENMP_FOUND)
    SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
  ENDIF()
ENDIF()

# - Find FFTW
# Find the native FFTW includes and library
#
#  FFTW_INCLUDES    - where to find fftw3.h
#  FFTW_LIBRARIES   - List of libraries when using FFTW.
#  FFTW_FOUND       - True if FFTW found.

#if (FFTW_INCLUDES)
  # Already in cache, be silent
#  set (FFTW_FIND_QUIETLY TRUE)
#endif (FFTW_INCLUDES)

#find_path (FFTW_INCLUDES fftw3.h)

#find_library (FFTW_LIBRARIES NAMES fftw3)

# handle the QUIETLY and REQUIRED arguments and set FFTW_FOUND to TRUE if
# all listed variables are TRUE
#include (FindPackageHandleStandardArgs)
#find_package_handle_standard_args (FFTW DEFAULT_MSG FFTW_LIBRARIES FFTW_INCLUDES)

#if(FFTW_FOUND)
#  add_definitions(-DFFTW)
#  target_link_libraries(ACFcalculator ${FFTW_LIBRARIES})
#endif ()

#SET(CMAKE_FIND_LIBRARY_SUFFIXES ".a")
#SET(BUILD_SHARED_LIBRARIES OFF)
#SET(CMAKE_EXE_LINKER_FLAGS "-static")

add_executable(traj2diffusion traj2diffusion.cpp)
target_link_libraries(traj2diffusion ${Boost_LIBRARIES})
