#
# ref: http://stackoverflow.com/questions/7724569/debug-vs-release-in-cmake
#
# moved
cmake_minimum_required(VERSION 2.8.9)

project(creek)


# SET(BOOST_ROOT /usr/local/boost_1_61_0/)
# find_package(Boost
#   1.61 EXACT
#   REQUIRED
#   COMPONENTS date_time system thread program_options
# )
include_directories(
  ./
  ../CTPL
  ../Kaskade
#  ${Boost_INCLUDE_DIRS}
  # ../tbb/include
  # /ssd/local/include
  )

# set(BOOST_INCLUDE_DIRS /usr/include)
# set(BOOST_LIBRARY_DIRS /usr/lib/x86_64-linux-gnu/)

# if(BOOST_FOUND)
#   message("BOOST_VERSION", ${Boost_VERSION})
#   message("BOOST_INCLUDE_DIRS", ${Boost_INCLUDE_DIRS})
# endif()

# in case of use on arm cross toolchain
if(CMAKE_ENV STREQUAL "arm")
  # add required libraries
  # all libraries are in /shared/hikey
  include_directories(
    /shared/hikey/include
    )
  link_directories(
    /shared/hikey/lib
    /shared/hikey/lib/aarch64-linux-gnu
    )
  SET(CMAKE_SYSTEM_NAME Linux)
  SET(CMAKE_SYSTEM_PROCESSOR aarch64)

  SET(COMPILER_PATH "/shared/hikey/toolchains/aarch64/bin")
  SET(TARGET_ABI "linux-gnu")
  # specify the cross compiler
  SET(CMAKE_C_COMPILER "${COMPILER_PATH}/aarch64-${TARGET_ABI}-gcc")
  SET(CMAKE_CXX_COMPILER "${COMPILER_PATH}/aarch64-${TARGET_ABI}-g++")

  SET(CMAKE_FIND_ROOT_PATH)
  # search for programs in the build host directories
  SET(CMAKE_FIND_ROOT_PATH_MODE_PROGRAM NEVER)
  # for libraries and headers in the target directories
  SET(CMAKE_FIND_ROOT_PATH_MODE_LIBRARY ONLY)
  SET(CMAKE_FIND_ROOT_PATH_MODE_INCLUDE ONLY)
  SET(CMAKE_FIND_ROOT_PATH_MODE_PACKAGE ONLY)

  SET(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -Wno-deprecated-declarations" ) # remove warning of auto_ptr
  #elseif(CMAKE_ENV STREQUAL "x86")
else()    # default: x86
  #SET( CMAKE_CXX_COMPILER "clang++" )  # has some issues with template template class...
  SET( CMAKE_CXX_COMPILER "g++" )
  #SET( CMAKE_CXX_COMPILER "icc" )
endif(CMAKE_ENV STREQUAL "arm")

# CLion needs to know all the headers.
# produced by ``ls -1"
SET(HEADERS
  AppliedPTransform.h
  Source/Bounded.h
  Source/BoundedCol2Evaulator.h
  Source/BoundedEvaulator.h
  BoundedWindow.h
  Create.h
  DoFn.h
  core/DualInputTransformEvaluator.h
  core/EvaluationBundleContext.h
  core/ExecutionContext.h
  core/Executor.h
  Win/FixedWindowInto.h
  Win/FixedWindowIntoEvaluator.h
  Win/GlobalWindow.h
  Mapper/GrepMapper.h
  Mapper/GrepMapperEvaluator.h
  Join/Join.h
  Join/JoinEvaluator.h
  Join/JoinEvaluator1.h
  LookAheadPipelineVisitor.h
  Mapper/Mapper.h
  MinimalContext.h
  MinimalProcessContext.h
  Source/NetworkLatencySource.h
  Source/NetworkLatencySourceEvaluator.h
  NopEvaluator.h
  OffsetRangeTracker.h
  PaneInfo.h
  ParDo/ParDo.h
  ParDo/ParDoCol2.h
  ParDo/ParDoEvaluator.h
  ParDo/ParDoSimpler.h
  ParDo/ParDoSimplerEvaluator.h
  ParDo/ParUpdateCol2.h
  ParDo/ParUpdateCol2Evaluator.h
  PerBundleContext.h
  core/Pipeline.h
  core/PipelineRunner.h
  core/PipelineVisitor.h
  Sink/RecordBitmapBundleSinkEvaluator.h
  Sink/RecordBitmapBundleSinkEvaluator1.h
  Sink/RecordBundleSinkEvaluator.h
  Sink/RecordBundleSinkEvaluatorJD.h
  Select/Select.h
  Win/SessionWindowInto.h
  Win/SessionWindowIntoEvaluator.h
  Mapper/SimpleMapper.h
  Mapper/SimpleMapperEvaluator.h
  Select/SimpleSelect.h
  Select/SimpleSelectEvaluator.h
  core/SingleInputTransformEvaluator.h
  Sink/Sink.h
  core/StatefulTransform.h
  StringUtils.h
  core/TransformEvaluator.h
  core/TransformTreeNode.h
  core/Transforms-multi.h
  core/Transforms.h
  Source/Unbounded-tweets.h
  Source/Unbounded.h
  Source/UnboundedInMemEvaluator.h
  Source/UnboundedInMemEvaluatorBase.h
  Source/UnboundedInMemEvaluator_2out.h
  Source/UnboundedInMemEvaluator_Join.h
  Source/UnboundedPseudoEvaluator.h
  Source/Unbounded_Join.h
  ValueContainer.h
  Values.h
  Win/WinGBK.h
  Win/WinGBKEvaluator.h
  Mapper/WinMapper.h
  Mapper/WinMapperEvaluator.h
  WinSum/WinSum-dispatch-eval.h
  WinSum/WinSumBase.h
  WinSum/WinSumEval.h
  WinSum/WinSum_addlong.h
  WinSum/WinSum_mergeset.h
  WinSum/WinSum_mergevector.h
  WinSum/WinSum_addlong_tvpair.h
  WinSum/WinSumTvEval.h
  WindowFn.h
  WinKeyFrag.h
  WinKeyReducer/WinKeyReducer.h
  WinKeyReducer/WinKeyReducerEval.h
  #WinKeyReducer/WinKeyReducerEval-MT.h
  #WinKeyReducer/WinKeyReducerEvalRecordBundle.h
  WinKeyReducer/WinKeyReducer-wc-common.h
  WinKeyReducer/WinKeyReducerEval-wc-common.h
  Mapper/WindowedGrepMapper.h
  Mapper/WindowedGrepMapperEvaluator.h
  WinSum/WindowedSumEvaluator.h
  WindowingStrategy.h
  Sink/WindowsBundleSinkEvaluator.h
  WindowsHolder.h
  Mapper/WordCountMapper.h
  Mapper/WordCountMapperEvaluator.h
  config.h
  creek-map.h
  creek-types.h
  linux-sizes.h
  log.h
  measure.h
  )

# common source shared
SET(SRC_WINSUM
  WinSum/WinSumBase.cpp
  WinSum/WinSum_addlong.cpp
  WinSum/WinSum_mergeset.cpp
  WinSum/WinSum_mergevector.cpp
  WinSum/WinSum_addlong_tvpair.cpp
  )

SET(SRC_WKR
  WinKeyReducer/WinKeyReducer.cpp
  WinKeyReducer/WinKeyReducer-wc-winbundle.cpp
  WinKeyReducer/WinKeyReducer-wc-recordbundle.cpp
  WinKeyReducer/WinKeyReducer-netmon.cpp
  WinKeyReducer/WinKeyReducerEval-wc-winbundle.cpp
  WinKeyReducer/WinKeyReducerEval-wc-recordbundle.cpp
  WinKeyReducer/WinKeyReducerEval-netmon.cpp
  )

SET(SOURCES
  ${HEADERS}
  core/Pipeline.cpp BoundedWindow.cpp Values.cpp core/Transform.cpp
  Source/Unbounded.cpp
  Mapper/WordCountMapper.cpp
  Win/WinGBK.cpp
  Mapper/GrepMapper.cpp
  Mapper/WindowedGrepMapper.cpp
  Sink/Sink.cpp
  Sink/WindowsBundleSink.cpp
  LookAheadPipelineVisitor.cpp
  Win/FixedWindowInto.cpp
  measure.c
  core/EvaluationBundleContext.cpp
  Mapper/SimpleMapper.cpp
  Join/Join.cpp
  Source/NetworkLatencySource.cpp
  #../Executor.cpp
  #../Kaskade/utilities/kalloc.cpp         # for numa allocator
  #../Kaskade/utilities/threading.cpp      # for numa allocator
  test/test-common.cpp
  Source/Unbounded_Join.cpp
  ${SRC_WINSUM}
  ${SRC_WKR}
  )

# -- apply to all configurations
#SET( CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -std=c++11 -g -fcilkplus -Wall" )
#SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -include config.h -fdiagnostics-color=auto -D_GLIBCXX_USE_CXX11_ABI=0")
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -include config.h -fdiagnostics-color=auto")

SET(CMAKE_CXX_FLAGS_DEBUG
  "${CMAKE_CXX_FLAGS_DEBUG} -std=c++11 -g -O0 -Wall -DHAVE_NUMA -Wfatal-errors")

#SET( CMAKE_CXX_FLAGS_RELEASE  "${CMAKE_CXX_FLAGS_RELEASE} -std=c++11 -g -fcilkplus -Wall -O2" )
SET(CMAKE_CXX_FLAGS_RELEASE "-std=c++11 -Wall -O3 -DHAVE_NUMA")

SET(CMAKE_CXX_FLAGS_MINSIZEREL
  "${CMAKE_CXX_FLAGS} -std=c++11 -g -fcilkplus -Wall")

if (CMAKE_BUILD_TYPE STREQUAL "Debug")
  add_definitions(-DCONFIG_KAGE_GLOBAL_DEBUG_LEVEL=30)
  #	add_definitions(-D_GLIBCXX_DEBUG)	   # will fail compilation in folly
endif ()

if (CMAKE_BUILD_TYPE STREQUAL "Release")
  add_definitions(-DCONFIG_KAGE_GLOBAL_DEBUG_LEVEL=50)
  add_definitions(-DNDEBUG)
endif ()

if (CMAKE_BUILD_TYPE STREQUAL "MinSizeRel")
  add_definitions(-DCONFIG_KAGE_GLOBAL_DEBUG_LEVEL=20)
endif ()

# SET( CMAKE_EXE_LINKER_FLAGS  "${CMAKE_EXE_LINKER_FLAGS} -lboost_date_time" )
SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -pthread")
#SET( CMAKE_EXE_LINKER_FLAGS  "${CMAKE_EXE_LINKER_FLAGS}" )
#SET( CMAKE_EXE_LINKER_FLAGS  "${CMAKE_EXE_LINKER_FLAGS} -pthread -L/home/xzl/tmp/boost_1_61_0/stage/lib -static" )

#add_executable(CreateTest.bin CreateTest.cpp ${SOURCES})

add_library(creek STATIC
  ${SOURCES}
  )

add_executable(test-grep.bin
  test/test-grep.cpp
  Sink/RecordBundleSink.cpp
  ${SOURCES}
  )

add_executable(test-join.bin test/test-join.cpp ${SOURCES})

add_executable(
  test-wingrep.bin test/test-wingrep.cpp
  Sink/RecordBundleSink.cpp
  ${SOURCES}
  )
add_executable(test-distinct.bin
  test/test-distinct.cpp
  ${SOURCES}
  Sink/RecordBundleSink.cpp
  )

add_executable(test-join-2.bin test/test-join-2.cpp ${SOURCES})

#add_executable(networklatency.bin
#        test/networklatency.cpp
#        core/Pipeline.cpp BoundedWindow.cpp Values.cpp core/Transform.cpp
#        Source/Unbounded.cpp
#        Win/WinGBK.cpp
#        #WindowKeyedReducer.cpp
#        WinKeyReducer/WinKeyReducer-netmon.cpp
#        #Sink.cpp
#        Sink/WindowsBundleSink.cpp
#        LookAheadPipelineVisitor.cpp
#        Win/FixedWindowInto.cpp
#        measure.c
#        core/EvaluationBundleContext.cpp
#        Source/NetworkLatencySource.cpp
#        #Kaskade/utilities/kalloc.cpp         # for numa allocator
#        #Kaskade/utilities/threading.cpp      # for numa allocator
#        test/test-common.cpp
#        WinKeyReducer/WinKeyReducerEval-netmon.cpp
#        )

add_executable(networklatency.bin
  test/networklatency.cpp
  )

add_executable(test-wc.bin
  test/test-wc.cpp
  )



add_executable(test-yahoo.bin
  test/test-yahoo.cpp
  Yahoo/YahooRecord.h
  Yahoo/YahooParser.h
  Yahoo/YahooParser.cpp
  Yahoo/YahooParserEvaluator.h
  Yahoo/YahooFilter.h
  Yahoo/YahooFilter.cpp
  Yahoo/YahooFilterEvaluator.h
  Yahoo/YahooAggregation.h
  Yahoo/YahooAggregation.cpp
  Yahoo/YahooAggregationEvaluator.h
  WinKeyReducer/WinKeyReducer-yahoo-common.h
  WinKeyReducer/WinKeyReducer-yahoo-winbundle.cpp
  # WinKeyReducer/WinKeyReducer-yahoo-recordbundle.cpp
  WinKeyReducer/WinKeyReducerEval-yahoo-common.h
  WinKeyReducer/WinKeyReducerEval-yahoo-winbundle.cpp
  # WinKeyReducer/WinKeyReducerEval-yahoo-recordbundle.cpp
  )

add_executable(test-yahoo-pc.bin
  test/test-yahoo.cpp
  Yahoo/YahooRecord.h
  Yahoo/YahooParser.h
  Yahoo/YahooParser.cpp
  Yahoo/YahooParserEvaluator.h
  Yahoo/YahooFilter.h
  Yahoo/YahooFilter.cpp
  Yahoo/YahooFilterEvaluator.h
  Yahoo/YahooAggregation.h
  Yahoo/YahooAggregation.cpp
  Yahoo/YahooAggregationEvaluator.h
  WinKeyReducer/WinKeyReducer-yahoo-common.h
  WinKeyReducer/WinKeyReducer-yahoo-winbundle.cpp
  # WinKeyReducer/WinKeyReducer-yahoo-recordbundle.cpp
  WinKeyReducer/WinKeyReducerEval-yahoo-common.h
  WinKeyReducer/WinKeyReducerEval-yahoo-winbundle.cpp
  papi/papi_sampler.hpp
  papi/papi_sampler.cpp
  # WinKeyReducer/WinKeyReducerEval-yahoo-recordbundle.cpp
  )



add_executable(test-yahoo-countbased.bin
  test/test-yahoo-countbased.cpp
  Yahoo/YahooRecord.h
  Yahoo/YahooParser.h
  Yahoo/YahooParser.cpp
  Yahoo/YahooParserEvaluator.h
  Yahoo/YahooFilter.h
  Yahoo/YahooFilter.cpp
  Yahoo/YahooFilterEvaluator.h
  Yahoo/YahooAggregation.h
  Yahoo/YahooAggregation.cpp
  Yahoo/YahooAggregationEvaluator.h
  Win/CountBasedTumblingWindow.h
  Win/CountBasedTumblingWindow.cpp
  Win/CountBasedTumblingWindowEvaluator.h

  WinKeyReducer/WinKeyReducer-yahoo-common.h
  WinKeyReducer/WinKeyReducer-yahoo-winbundle.cpp
  # WinKeyReducer/WinKeyReducer-yahoo-recordbundle.cpp
  WinKeyReducer/WinKeyReducerEval-yahoo-common.h
  WinKeyReducer/WinKeyReducerEval-yahoo-winbundle.cpp
  # WinKeyReducer/WinKeyReducerEval-yahoo-recordbundle.cpp
  )

add_executable(test-yahoo-countbased-pc.bin
  test/test-yahoo-countbased.cpp
  Yahoo/YahooRecord.h
  Yahoo/YahooParser.h
  Yahoo/YahooParser.cpp
  Yahoo/YahooParserEvaluator.h
  Yahoo/YahooFilter.h
  Yahoo/YahooFilter.cpp
  Yahoo/YahooFilterEvaluator.h
  Yahoo/YahooAggregation.h
  Yahoo/YahooAggregation.cpp
  Yahoo/YahooAggregationEvaluator.h
  Win/CountBasedTumblingWindow.h
  Win/CountBasedTumblingWindow.cpp
  Win/CountBasedTumblingWindowEvaluator.h

  WinKeyReducer/WinKeyReducer-yahoo-common.h
  WinKeyReducer/WinKeyReducer-yahoo-winbundle.cpp
  # WinKeyReducer/WinKeyReducer-yahoo-recordbundle.cpp
  WinKeyReducer/WinKeyReducerEval-yahoo-common.h
  WinKeyReducer/WinKeyReducerEval-yahoo-winbundle.cpp
  # WinKeyReducer/WinKeyReducerEval-yahoo-recordbundle.cpp
  papi/papi_sampler.hpp
  papi/papi_sampler.cpp
  )


add_executable(test-yahoo-lsds.bin
  test/test-yahoo-lsds.cpp
  Source/YahooBenchmarkSource.cpp
  Source/YahooBenchmarkSource.h
  Source/YahooBenchmarkSourceEvaluator.h
  Mapper/YahooMapper.cpp
  Mapper/YahooMapper.h
  Mapper/YahooMapperEvaluator.h
  Sink/RecordBundleSink.cpp
  Sink/WindowsBundleSink.cpp
  )


add_executable(test-yahoo-lsds1.bin
  test/test-yahoo-lsds1.cpp
  Source/YahooBenchmarkSource.cpp
  Source/YahooBenchmarkSource.h
  Source/YahooBenchmarkSourceEvaluator.h
  Mapper/YahooMapper.cpp
  Mapper/YahooMapper.h
  Mapper/YahooMapperEvaluator.h
  Sink/RecordBundleSink.cpp
  Sink/WindowsBundleSink.cpp
  )

add_executable(test-yahoo-lsds2.bin
  test/test-yahoo-lsds2.cpp
  # Source/YahooBenchmarkSource.cpp
  # Source/YahooBenchmarkSource.h
  # Source/YahooBenchmarkSourceEvaluator.h
  Yahoo/YahooRecord.h
  Yahoo/YahooParser.h
  Yahoo/YahooParser.cpp
  Yahoo/YahooParserEvaluator.h

  Mapper/YahooMapper.cpp
  Mapper/YahooMapper.h
  Mapper/YahooMapperEvaluator.h
  Sink/RecordBundleSink.cpp
  Sink/WindowsBundleSink.cpp
  )

add_executable(test-nyt.bin
  test/test-nyt.cpp
  NewYorkTaxi/NYTRecord.h
  NewYorkTaxi/NYTParser.h
  NewYorkTaxi/NYTParser.cpp
  NewYorkTaxi/NYTParserEvaluator.h
  NewYorkTaxi/NYTQuery.h
  NewYorkTaxi/NYTQuery.cpp
  NewYorkTaxi/NYTQueryEvaluator.h
  NewYorkTaxi/NYTWinKeyReducerCommon.h
  NewYorkTaxi/NYTWinKeyReducerEvalCommon.h
  NewYorkTaxi/NYTWinKeyReducerEvalWinBundle.cpp
  NewYorkTaxi/NYTWinKeyReducerWinBundle.cpp

  )

add_executable(test-nyt-pc.bin
  test/test-nyt.cpp
  NewYorkTaxi/NYTRecord.h
  NewYorkTaxi/NYTParser.h
  NewYorkTaxi/NYTParser.cpp
  NewYorkTaxi/NYTParserEvaluator.h
  NewYorkTaxi/NYTQuery.h
  NewYorkTaxi/NYTQuery.cpp
  NewYorkTaxi/NYTQueryEvaluator.h
  NewYorkTaxi/NYTWinKeyReducerCommon.h
  NewYorkTaxi/NYTWinKeyReducerEvalCommon.h
  NewYorkTaxi/NYTWinKeyReducerEvalWinBundle.cpp
  NewYorkTaxi/NYTWinKeyReducerWinBundle.cpp
  papi/papi_sampler.cpp
  papi/papi_sampler.hpp
  )


add_executable(test-window.bin
  test/test-window.cpp
  NewYorkTaxi/NYTRecord.h
  NewYorkTaxi/NYTParser.h
  NewYorkTaxi/NYTParser.cpp
  NewYorkTaxi/NYTParserEvaluator.h
  NewYorkTaxi/NYTQuery.h
  NewYorkTaxi/NYTQuery.cpp
  NewYorkTaxi/NYTQueryEvaluator.h
  NewYorkTaxi/NYTWinKeyReducerCommon.h
  NewYorkTaxi/NYTWinKeyReducerEvalCommon.h
  NewYorkTaxi/NYTWinKeyReducerEvalWinBundle.cpp
  NewYorkTaxi/NYTWinKeyReducerWinBundle.cpp
  NewYorkTaxi/NYTWinKeyReducerEvalWinBundle.cpp
  NewYorkTaxi/NYTWinKeyReducerWinBundle.cpp

  NewYorkTaxi/NYTWindowedSumEvaluator.h
  NewYorkTaxi/NYTWindowedSum.h
  NewYorkTaxi/NYTWindowedSum.cpp
  Win/CountBasedTumblingWindow.h
  Win/CountBasedTumblingWindow.cpp
  Win/CountBasedTumblingWindowEvaluator.h
  )

add_executable(test-nyt-countbased.bin
  test/test-nyt-countbased.cpp
  NewYorkTaxi/NYTRecord.h
  NewYorkTaxi/NYTParser.h
  NewYorkTaxi/NYTParser.cpp
  NewYorkTaxi/NYTParserEvaluator.h
  NewYorkTaxi/NYTQuery.h
  NewYorkTaxi/NYTQuery.cpp
  NewYorkTaxi/NYTQueryEvaluator.h
  NewYorkTaxi/NYTWinKeyReducerCommon.h
  NewYorkTaxi/NYTWinKeyReducerEvalCommon.h
  NewYorkTaxi/NYTWinKeyReducerEvalWinBundle.cpp
  NewYorkTaxi/NYTWinKeyReducerWinBundle.cpp
  NewYorkTaxi/NYTWinKeyReducerEvalWinBundle.cpp
  NewYorkTaxi/NYTWinKeyReducerWinBundle.cpp

  NewYorkTaxi/NYTWindowedSumEvaluator.h
  NewYorkTaxi/NYTWindowedSum.h
  NewYorkTaxi/NYTWindowedSum.cpp
  Win/CountBasedTumblingWindow.h
  Win/CountBasedTumblingWindow.cpp
  Win/CountBasedTumblingWindowEvaluator.h
  )

add_executable(test-nyt-countbased-pc.bin
  test/test-nyt-countbased.cpp
  NewYorkTaxi/NYTRecord.h
  NewYorkTaxi/NYTParser.h
  NewYorkTaxi/NYTParser.cpp
  NewYorkTaxi/NYTParserEvaluator.h
  NewYorkTaxi/NYTQuery.h
  NewYorkTaxi/NYTQuery.cpp
  NewYorkTaxi/NYTQueryEvaluator.h
  NewYorkTaxi/NYTWinKeyReducerCommon.h
  NewYorkTaxi/NYTWinKeyReducerEvalCommon.h
  NewYorkTaxi/NYTWinKeyReducerEvalWinBundle.cpp
  NewYorkTaxi/NYTWinKeyReducerWinBundle.cpp
  NewYorkTaxi/NYTWinKeyReducerEvalWinBundle.cpp
  NewYorkTaxi/NYTWinKeyReducerWinBundle.cpp

  NewYorkTaxi/NYTWindowedSumEvaluator.h
  NewYorkTaxi/NYTWindowedSum.h
  NewYorkTaxi/NYTWindowedSum.cpp
  Win/CountBasedTumblingWindow.h
  Win/CountBasedTumblingWindow.cpp
  Win/CountBasedTumblingWindowEvaluator.h
  papi/papi_sampler.hpp
  papi/papi_sampler.cpp
  )

target_compile_definitions(creek
  PUBLIC LSDS_STREAMBOX=1
  # PUBLIC MEASURE_40M_TUPLES=0
  )

target_compile_definitions(test-yahoo-lsds.bin
  PUBLIC LSDS_STREAMBOX=1
  )

target_compile_definitions(test-yahoo-lsds1.bin
  PUBLIC LSDS_STREAMBOX=1
  )

target_compile_definitions(test-yahoo-lsds2.bin
  PUBLIC LSDS_STREAMBOX=1
  )

target_compile_definitions(test-yahoo-pc.bin
  PUBLIC MEASURE_40M_TUPLES=1
  )

#add_executable(test-wc.bin
#        test/test-wc.cpp
#        core/Pipeline.cpp BoundedWindow.cpp Values.cpp core/Transform.cpp
#        Source/Unbounded.cpp
#        Mapper/WordCountMapper.cpp
#        Win/WinGBK.cpp
#        #WindowKeyedReducer.cpp
#        WinKeyReducer/WinKeyReducer-wc.cpp
#        ${SRC_WINSUM}
#        WinKeyReducer/WinKeyReducer.cpp
#        Mapper/GrepMapper.cpp
#        Mapper/WindowedGrepMapper.cpp
#        #Sink.cpp
#        Sink/WindowsBundleSink.cpp
#        LookAheadPipelineVisitor.cpp
#        Win/FixedWindowInto.cpp
#        measure.c
#        core/EvaluationBundleContext.cpp
#        #SimpleMapper.cpp
#        #Join.cpp
#        #NetworkLatencySource.cpp
#        #Executor.cpp
#        #Kaskade/utilities/kalloc.cpp         # for numa allocator
#        #Kaskade/utilities/threading.cpp      # for numa allocator
#        test/test-common.cpp
#        WinKeyReducer/WinKeyReducerEval-wc.cpp
#        #WindowKeyedReducerEval-netmon.cpp
#        )

add_executable(test-tweet.bin
  test/test-tweet.cpp
  core/Pipeline.cpp BoundedWindow.cpp Values.cpp core/Transform.cpp
  Source/Unbounded.cpp
  Mapper/WordCountMapper.cpp
  Win/WinGBK.cpp
  Mapper/WinMapper.cpp
  #WindowKeyedReducer.cpp
  WinKeyReducer/WinKeyReducer-wc-winbundle.cpp
  WinKeyReducer/WinKeyReducer-wc-recordbundle.cpp
  WinKeyReducer/WinKeyReducerEval-wc-winbundle.cpp
  WinKeyReducer/WinKeyReducerEval-wc-recordbundle.cpp
  #WindowKeyedReducer-netmon.cpp
  #WinSum/WindowedSumEvaluator1-tvpair.cpp
  ${SRC_WINSUM}
  #WindowKeyedReducer.cpp
  #GrepMapper.cpp
  #WindowedGrepMapper.cpp
  Sink/Sink.cpp
  Sink/WindowsBundleSink.cpp
  Sink/RecordBundleSink.cpp
  LookAheadPipelineVisitor.cpp
  #FixedWindowInto.cpp
  measure.c
  core/EvaluationBundleContext.cpp
  #SimpleMapper.cpp
  Join/Join.cpp
  #NetworkLatencySource.cpp
  #Executor.cpp
  #Kaskade/utilities/kalloc.cpp         # for numa allocator
  #Kaskade/utilities/threading.cpp      # for numa allocator
  test/test-common.cpp
  #WindowKeyedReducerEval-wc.cpp
  #WindowKeyedReducerEval-netmon.cpp
  )

set_target_properties(test-tweet.bin
  #PROPERTIES COMPILE_FLAGS "-DWORKAROUND_JOIN_JDD=1 -DWORKAROUND_WINKEYREDUCER_RECORDBUNDLE=1"
  PROPERTIES COMPILE_FLAGS "-DWORKAROUND_JOIN_JDD=1"
  )

#add_executable(test-WordCountMapper.bin test-WordCountMapper.cpp ${SOURCES})

# -- uncomment to turn on the target
#add_executable(test-kaskade test-kaskade.cpp
#Kaskade/utilities/kalloc.cpp
#Kaskade/utilities/threading.cpp
#)

# xzl: tbb will be in the default lib path
#FIND_LIBRARY(tbb NAMES libtbb.so)
FIND_LIBRARY(tbb NAMES libtbb.so PATHS tbb/lib/intel64/gcc4.7)
FIND_LIBRARY(papi NAMES libpapi.so)
#FIND_LIBRARY(tbb_malloc NAMES libtbbmalloc.so PATHS tbb/lib/intel64/gcc4.7)

# x86 libs
if(NOT CMAKE_ENV STREQUAL "arm")
  FIND_LIBRARY(folly NAMES libfolly.a PATHS /usr/local/lib)
  FIND_LIBRARY(re2 NAMES libre2.a PATHS /usr/local/lib)
  FIND_LIBRARY(glog NAMES libglog.a PATHS /usr/local/lib)
  FIND_LIBRARY(re2 NAMES libre2.a PATHS /usr/local/lib)
  FIND_LIBRARY(jemalloc NAMES libjemalloc.a PATHS /usr/local/lib)
  #FIND_LIBRARY(cityhash NAMES libcityhash.a PATHS /ssd/local/lib)
endif()

#TARGET_LINK_LIBRARIES(CreateTest.bin boost_date_time boost_system boost_thread numa ${tbb})
#TARGET_LINK_LIBRARIES(test-join.bin boost_date_time boost_system boost_thread numa ${tbb})

TARGET_LINK_LIBRARIES(test-join-2.bin
  boost_date_time
  boost_system
  boost_thread
  boost_program_options
  numa ${tbb}
  ${folly}
  ${re2}
  pcre
  )


TARGET_LINK_LIBRARIES(networklatency.bin
  creek
  boost_date_time
  boost_system
  boost_thread
  boost_program_options
  numa ${tbb}
  ${folly} ${glog} double-conversion  # needed by atomic hash map
  ${re2}
  pcre
  )

TARGET_LINK_LIBRARIES(test-join.bin
  boost_date_time
  boost_system
  boost_thread
  boost_program_options
  numa ${tbb}
  ${folly}
  ${re2}
  pcre
  )

TARGET_LINK_LIBRARIES(test-wc.bin
  creek
  boost_date_time
  boost_system
  boost_thread
  boost_program_options
  numa
  ${tbb}
  #jemalloc -- do we have this for ARM?
  #  ${tbb_malloc}
  ${folly}
  ${re2}
  pcre
  )

TARGET_LINK_LIBRARIES(test-yahoo.bin
  creek
  boost_date_time
  boost_system
  boost_thread
  boost_program_options
  numa
  ${tbb}
  #jemalloc -- do we have this for ARM?
  #  ${tbb_malloc}
  ${folly}
  pcre
  )



TARGET_LINK_LIBRARIES(test-yahoo-countbased.bin
  creek
  boost_date_time
  boost_system
  boost_thread
  boost_program_options
  numa
  ${tbb}
  #jemalloc -- do we have this for ARM?
  #  ${tbb_malloc}
  ${folly}
  pcre
  )


TARGET_LINK_LIBRARIES(test-yahoo-lsds.bin
  creek
  boost_date_time
  boost_system
  boost_thread
  boost_program_options
  boost_regex
  boost_iostreams
  numa
  ${tbb}
  #jemalloc -- do we have this for ARM?
  #  ${tbb_malloc}
  ${folly}
  ${re2}
  pcre
  )

TARGET_LINK_LIBRARIES(test-yahoo-lsds1.bin
  creek
  boost_date_time
  boost_system
  boost_thread
  boost_program_options
  boost_regex
  boost_iostreams
  numa
  ${tbb}
  #jemalloc -- do we have this for ARM?
  #  ${tbb_malloc}
  ${folly}
  ${re2}
  pcre
  )


TARGET_LINK_LIBRARIES(test-yahoo-lsds2.bin
  creek
  boost_date_time
  boost_system
  boost_thread
  boost_program_options
  boost_regex
  boost_iostreams
  numa
  ${tbb}
  #jemalloc -- do we have this for ARM?
  #  ${tbb_malloc}
  ${folly}
  ${re2}
  pcre
  )

# target_compile_definitions(test-yahoo-lsds.bin
#   PUBLIC LSDS_STREAMBOX=1)
TARGET_LINK_LIBRARIES(test-tweet.bin
  boost_date_time
  boost_system
  boost_thread
  boost_program_options
  numa ${tbb}
  ${jemalloc}
  ${folly}
  ${re2}
  pcre
  )

TARGET_LINK_LIBRARIES(test-grep.bin
  boost_date_time
  boost_system
  boost_thread
  boost_program_options
  numa ${tbb}
  ${folly}
  ${re2}
  pcre
  )

TARGET_LINK_LIBRARIES(test-wingrep.bin
  boost_date_time
  boost_system
  boost_thread
  boost_program_options
  numa ${tbb}
  ${folly}
  ${re2}
  pcre
  )

TARGET_LINK_LIBRARIES(test-distinct.bin
  boost_date_time
  boost_system
  boost_thread
  boost_program_options
  numa ${tbb}
  ${folly} ${glog} double-conversion  # needed by atomic hash map
  ${cityhash} # needed by libcuckoo
  ${re2}
  pcre
  )

#TARGET_LINK_LIBRARIES(test-WordCountMapper.bin
#	boost_date_time
#	boost_system
#	boost_thread
#	numa ${tbb}
#	re2
#	pcre
#)
#TARGET_LINK_LIBRARIES(test-kaskade boost_date_time boost_system boost_thread numa)

# can't name test
#add_custom_target(all
#	DEPENDS test-wc.bin test-grep.bin test-distinct.bin test-wingrep.bin
#	)

add_executable(test-linear-road.bin
  test/test-linear-road.cpp
  LinearRoad/LinearRoadParser.h
  LinearRoad/LinearRoadParser.cpp
  LinearRoad/LinearRoadParserEvaluator.h
  LinearRoad/LinearRoadRecord.h
  LinearRoad/LinearRoadFilter.h
  LinearRoad/LinearRoadFilter.cpp
  LinearRoad/LinearRoadFilterEvaluator.h
  LinearRoad/LinearRoadTollNotification.h
  LinearRoad/LinearRoadTollNotification.cpp
  LinearRoad/LinearRoadTollNotificationEvaluator.h

  WinKeyReducer/WinKeyReducer-yahoo-common.h
  WinKeyReducer/WinKeyReducer-yahoo-winbundle.cpp
  # WinKeyReducer/WinKeyReducer-yahoo-recordbundle.cpp
  WinKeyReducer/WinKeyReducerEval-yahoo-common.h
  WinKeyReducer/WinKeyReducerEval-yahoo-winbundle.cpp
  )

add_executable(test-nexmark.bin
  test/test-nexmark.cpp
  Nexmark/NexmarkRecord.hpp
  Nexmark/NexmarkParser.hpp
  Nexmark/NexmarkParser.cpp
  Nexmark/NexmarkParserEvaluator.hpp
  Nexmark/NexmarkFilter.hpp
  Nexmark/NexmarkFilter.cpp
  Nexmark/NexmarkFilterEvaluator.hpp
  Nexmark/NexmarkAggregation.hpp
  Nexmark/NexmarkAggregation.cpp
  Nexmark/NexmarkAggregationEvaluator.hpp

  )

TARGET_LINK_LIBRARIES(test-nexmark.bin
  creek
  boost_date_time
  boost_system
  boost_thread
  boost_program_options
  numa ${tbb}
  ${folly}
  ${re2}
  pcre
  )

TARGET_LINK_LIBRARIES(test-linear-road.bin
  creek
  boost_date_time
  boost_system
  boost_thread
  boost_program_options
  numa ${tbb}
  ${folly}
  ${re2}
  pcre
  )

TARGET_LINK_LIBRARIES(test-nyt.bin
  creek
  boost_date_time
  boost_system
  boost_thread
  boost_program_options
  numa
  ${tbb}
  #jemalloc -- do we have this for ARM?
  #  ${tbb_malloc}
  ${folly}
  pcre
  )


TARGET_LINK_LIBRARIES(test-window.bin
  creek
  boost_date_time
  boost_system
  boost_thread
  boost_program_options
  numa
  ${tbb}
  #jemalloc -- do we have this for ARM?
  #  ${tbb_malloc}
  ${folly}
  pcre
  )

TARGET_LINK_LIBRARIES(test-nyt-countbased.bin
  creek
  boost_date_time
  boost_system
  boost_thread
  boost_program_options
  numa
  ${tbb}
  #jemalloc -- do we have this for ARM?
  #  ${tbb_malloc}
  ${folly}
  pcre
  )


TARGET_LINK_LIBRARIES(test-yahoo-pc.bin
  creek
  boost_date_time
  boost_system
  boost_thread
  boost_program_options
  numa
  ${tbb}
  #jemalloc -- do we have this for ARM?
  #  ${tbb_malloc}
  ${folly}
  pcre
  papi
  )


TARGET_LINK_LIBRARIES(test-yahoo-countbased-pc.bin
  creek
  boost_date_time
  boost_system
  boost_thread
  boost_program_options
  numa
  ${tbb}
  #jemalloc -- do we have this for ARM?
  #  ${tbb_malloc}
  ${folly}
  pcre
  papi
  )

TARGET_LINK_LIBRARIES(test-nyt-pc.bin
  creek
  boost_date_time
  boost_system
  boost_thread
  boost_program_options
  numa
  ${tbb}
  #jemalloc -- do we have this for ARM?
  #  ${tbb_malloc}
  ${folly}
  pcre
  papi
  )

TARGET_LINK_LIBRARIES(test-nyt-countbased-pc.bin
  creek
  boost_date_time
  boost_system
  boost_thread
  boost_program_options
  numa
  ${tbb}
  #jemalloc -- do we have this for ARM?
  #  ${tbb_malloc}
  ${folly}
  pcre
  papi
  )
