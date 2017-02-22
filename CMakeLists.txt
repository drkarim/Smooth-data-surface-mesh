cmake_minimum_required(VERSION 2.8)
 
PROJECT(smoothscalars)

find_package(VTK REQUIRED)
include(${VTK_USE_FILE})


add_executable(smoothscalars SmoothScalars.cxx)
 
if(VTK_LIBRARIES)
  target_link_libraries(smoothscalars ${VTK_LIBRARIES})
  
endif()