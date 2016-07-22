cmake_minimum_required(VERSION 3.0)

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

