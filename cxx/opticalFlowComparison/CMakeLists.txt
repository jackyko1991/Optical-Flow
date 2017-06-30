cmake_minimum_required(VERSION 2.8)

project( OpticalFlowComparison )

find_package(ITK REQUIRED)
include(${ITK_USE_FILE})

find_package( OpenCV REQUIRED )

if(NOT MATLAB_FOUND)
  find_package(MATLAB REQUIRED)
  include_directories(${Matlab_INCLUDE_DIRS})
endif()

if (ITKVtkGlue_LOADED)
  find_package(VTK REQUIRED)
  include(${VTK_USE_FILE})
else()
  find_package(ItkVtkGlue REQUIRED)
  include(${ItkVtkGlue_USE_FILE})
  set(Glue ItkVtkGlue)
endif()

include_directories(${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR})
 
file(GLOB HEADER_FILES *.h)
file(GLOB CPP_FILES *.cpp)
file(GLOB CXX_FILES *.cxx)

add_executable( OpticalFlowComparison ${CXX_FILES} ${CPP_FILES} ${HEADER_FILES})
target_link_libraries( OpticalFlowComparison ${OpenCV_LIBS} ${Glue}  ${VTK_LIBRARIES} ${ITK_LIBRARIES} ${BuildMex_LIBRARIES})