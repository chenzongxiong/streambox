#
# ref: http://stackoverflow.com/questions/7724569/debug-vs-release-in-cmake
#
cmake_minimum_required(VERSION 2.8.9)

project(creek-micro)

include_directories(
../
../Kaskade
)

#SET( SOURCES
#../measure.c
#../Kaskade/utilities/kalloc.cpp
#../Kaskade/utilities/threading.cpp
#)

#SET( CMAKE_CXX_COMPILER "clang++" )
SET( CMAKE_CXX_COMPILER "g++-5" )
#SET( CMAKE_CXX_COMPILER "icc" )

# -- apply to all configurations
SET( CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -std=c++11 -Wall -g" )

#SET( CMAKE_CXX_FLAGS_DEBUG
#	"${CMAKE_CXX_FLAGS_DEBUG} -std=c++11 -g -O0 -Wall -DHAVE_NUMA" )

#SET( CMAKE_CXX_FLAGS_RELEASE  "${CMAKE_CXX_FLAGS_RELEASE} -std=c++11 -g -fcilkplus -Wall -O2" )
#SET( CMAKE_CXX_FLAGS_RELEASE  "-std=c++11 -g -fcilkplus -Wall -O2" )
#SET( CMAKE_CXX_FLAGS_RELEASE  "-std=c++11 -g -Wall -O2 -DHAVE_NUMA" )

#SET( CMAKE_CXX_FLAGS_MINSIZEREL
#			"${CMAKE_CXX_FLAGS} -std=c++11 -g -fcilkplus -Wall" )

if(CMAKE_BUILD_TYPE STREQUAL "Debug")
	add_definitions(-DCONFIG_KAGE_GLOBAL_DEBUG_LEVEL=30)
#	add_definitions(-DUSE_CILK)
	add_definitions(-D_GLIBCXX_DEBUG)
endif()

if(CMAKE_BUILD_TYPE STREQUAL "Release")
	add_definitions(-DCONFIG_KAGE_GLOBAL_DEBUG_LEVEL=50)
#	add_definitions(-DUSE_CILK)
endif()

if(CMAKE_BUILD_TYPE STREQUAL "MinSizeRel")
	add_definitions(-DCONFIG_KAGE_GLOBAL_DEBUG_LEVEL=20)
endif()

# SET( CMAKE_EXE_LINKER_FLAGS  "${CMAKE_EXE_LINKER_FLAGS} -lboost_date_time" )
SET( CMAKE_EXE_LINKER_FLAGS  "${CMAKE_EXE_LINKER_FLAGS} -pthread" )
#SET( CMAKE_EXE_LINKER_FLAGS  "${CMAKE_EXE_LINKER_FLAGS}" )

#add_executable(CreateTest.bin CreateTest.cpp ${SOURCES})

FIND_LIBRARY(tbb NAMES libtbb.so PATHS ../tbb/lib/intel64/gcc4.7)

# -- uncomment to turn on the target
add_executable(numa-tp.bin numa-tp.cpp
../Kaskade/utilities/kalloc.cpp
../Kaskade/utilities/threading.cpp
)

add_executable(mt-numa.bin mt-numa.cpp
../Kaskade/utilities/kalloc.cpp
../Kaskade/utilities/threading.cpp
)

TARGET_LINK_LIBRARIES(numa-tp.bin boost_thread boost_system numa)
TARGET_LINK_LIBRARIES(mt-numa.bin 
  boost_thread 
  boost_system 
  numa pcre
  ${tbb}
  )
