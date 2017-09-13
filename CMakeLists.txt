cmake_minimum_required (VERSION 2.8.11)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake" )

find_package(PCL 1.7 REQUIRED COMPONENTS common visualization io features registration )
find_package(GeographicLib 1.34 REQUIRED)
find_package(PCAP)
include_directories(${PCL_INCLUDE_DIRS} ${GeographicLib_INCLUDE_DIRS} ${PCAP_INCLUDE_DIR})
link_directories(${PCL_LIBRARY_DIRS})
add_definitions(${PCL_DEFINITIONS} ${GeographicLib_DEFINITIONS})


add_executable (pcdviewer SimplePCDViewer.cpp )
target_link_libraries ( pcdviewer
  ${PCL_COMMON_LIBRARIES}
  ${PCL_VISUALIZATION_LIBRARIES}
  ${PCL_IO_LIBRARIES}
  ${PCL_FEATURES_LIBRARIES}
  ${GeographicLib_LIBRARIES}
)

add_executable (pcapdump PcapDump.cpp MyHdlGrabber.cpp )
target_link_libraries ( pcapdump
  ${PCAP_LIBRARY}
  libpcap.a
  ${PCL_COMMON_LIBRARIES}
  ${PCL_IO_LIBRARIES}
  ${PCL_FEATURES_LIBRARIES}
  ${PCL_VISUALIZATION_LIBRARIES}
)

add_executable (cloud_reg_adapt CloudRegister.cpp)
SET_TARGET_PROPERTIES(cloud_reg_adapt PROPERTIES
  COMPILE_FLAGS "-DADAPT_SIMPLEX"
)

target_link_libraries ( cloud_reg_adapt
  ${PCL_COMMON_LIBRARIES}
  ${PCL_IO_LIBRARIES}
  ${PCL_FEATURES_LIBRARIES}
  ${PCL_VISUALIZATION_LIBRARIES}
  ${PCL_REGISTRATION_LIBRARIES}
  ${GeographicLib_LIBRARIES}
)

add_executable (lum LUM.cpp)
target_link_libraries ( lum
  ${PCL_COMMON_LIBRARIES}
  ${PCL_IO_LIBRARIES}
  ${PCL_FEATURES_LIBRARIES}
  ${PCL_VISUALIZATION_LIBRARIES}
  ${PCL_REGISTRATION_LIBRARIES}
  ${GeographicLib_LIBRARIES}
)
