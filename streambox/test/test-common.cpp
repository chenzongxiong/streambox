#include <stdio.h>
#include "log.h"
#include "config.h"

void print_config()
{
	printf("config:\n");

#ifdef USE_NUMA_TP
	xzl_show_define(USE_NUMA_TP);
#else
	xzl_show_undefine(USE_NUMA_TP);
#endif

#ifdef USE_TBB_DS
	xzl_show_define(USE_TBB_DS)
#else
	xzl_show_undefine(USE_TBB_DS)
#endif

#ifdef USE_FOLLY_VECTOR
	xzl_show_define(USE_FOLLY_VECTOR)
#else
	xzl_show_undefine(USE_FOLLY_VECTOR)
#endif

#ifdef USE_FOLLY_STRING
	xzl_show_define(USE_FOLLY_STRING)
#else
	xzl_show_undefine(USE_FOLLY_STRING)
#endif

#ifdef USE_FOLLY_HASHMAP
	xzl_show_define(USE_FOLLY_HASHMAP)
#else
	xzl_show_undefine(USE_FOLLY_HASHMAP)
#endif

#ifdef USE_CUCKOO_HASHMAP
	xzl_show_define(USE_CUCKOO_HASHMAP)
#else
	xzl_show_undefine(USE_CUCKOO_HASHMAP)
#endif

#ifdef USE_TBB_HASHMAP
	xzl_show_define(USE_TBB_HASHMAP)
#else
	xzl_show_undefine(USE_TBB_HASHMAP)
#endif

#ifdef USE_NUMA_ALLOC
	xzl_show_define(USE_NUMA_ALLOC)
#else
	xzl_show_undefine(USE_NUMA_ALLOC)
#endif

#ifdef INPUT_ALWAYS_ON_NODE0
	xzl_show_define(INPUT_ALWAYS_ON_NODE0)
#else
	xzl_show_undefine(INPUT_ALWAYS_ON_NODE0)
#endif

#ifdef MEASURE_LATENCY
	xzl_show_define(MEASURE_LATENCY)
#else
	xzl_show_undefine(MEASURE_LATENCY)
#endif

#ifdef CONFIG_SOURCE_THREADS
	xzl_show_define(CONFIG_SOURCE_THREADS)
#else
//	xzl_show_undefine(CONFIG_SOURCE_THREADS)
#error
#endif

#ifdef CONFIG_MIN_PERNODE_BUFFER_SIZE_MB
	xzl_show_define(CONFIG_MIN_PERNODE_BUFFER_SIZE_MB)
#else
//	xzl_show_undefine(CONFIG_MIN_PERNODE_BUFFER_SIZE_MB)
#error
#endif

#ifdef CONFIG_MIN_EPOCHS_PERNODE_BUFFER
	xzl_show_define(CONFIG_MIN_EPOCHS_PERNODE_BUFFER)
#else
//	xzl_show_undefine(X)
#error
#endif

#ifdef CONFIG_NETMON_HT_PARTITIONS
	xzl_show_define(CONFIG_NETMON_HT_PARTITIONS)
#else
#error
#endif

#ifdef CONFIG_JOIN_HT_PARTITIONS
	xzl_show_define(CONFIG_JOIN_HT_PARTITIONS)
#else
#error
#endif

	/* warn about workarounds */
#ifdef WORKAROUND_JOIN_JDD
	EE("warning: WORKAROUND_JOIN_JDD = 1. Sure? (any key to continue)\n");
	getchar();
#endif

#ifdef WORKAROUND_WINKEYREDUCER_RECORDBUNDLE
	EE("warning: WORKAROUND_WINKEYREDUCER_RECORDBUNDLE = 1. Sure? (any key to continue)\n");
	getchar();
#endif

	/* todo: add more here */

}

// Copyright Vladimir Prus 2002-2004.
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt
// or copy at http://www.boost.org/LICENSE_1_0.txt)

/* The simplest usage of the library.
 *
 * http://www.boost.org/doc/libs/1_61_0/libs/program_options/example/first.cpp
 * http://www.boost.org/doc/libs/1_61_0/libs/program_options/example/options_description.cpp
 */
#undef _GLIBCXX_DEBUG /* does not get along with program options lib */

#include <iostream>
#include <boost/program_options.hpp>
using namespace std;
namespace po = boost::program_options;

#include <thread>
#include "test-common.h"

void parse_options(int ac, char *av[], pipeline_config* config)
{
  po::variables_map vm;
  xzl_assert(config);

  try {

      po::options_description desc("Allowed options");
      desc.add_options()
          ("help", "produce help message")
          ("records", po::value<unsigned long>(), "records per wm interval")
          ("target_tput", po::value<unsigned long>(), "target throughput (krec/s)")
          ("record_size", po::value<unsigned long>(), "record (string_range) size (bytes)")
          ("window_size", po::value<unsigned long>(), "window_size (seconds)")
          ("window_count", po::value<unsigned long>(), "counts per window")
          ("campaigns", po::value<unsigned long>(), "campaigns")
          ("cores", po::value<unsigned long>(), "# cores for worker threads")
          ("input_file", po::value<vector<string>>(), "input file path") /* must be vector */
      ;

      po::store(po::parse_command_line(ac, av, desc), vm);
      po::notify(vm);

      if (vm.count("help")) {
          cout << desc << "\n";
          exit(1);
      }

      if (vm.count("records")) {
          config->records_per_interval = vm["records"].as<unsigned long>();
      }
      if (vm.count("target_tput")) {
          config->target_tput = vm["target_tput"].as<unsigned long>() * 1000;
      }
      if (vm.count("record_size")) {
          config->record_size = vm["record_size"].as<unsigned long>();
      }
      if (vm.count("window_size")) {
          config->window_size = vm["window_size"].as<unsigned long>();
      }
      if (vm.count("window_count")) {
          config->window_count = vm["window_count"].as<unsigned long>();
      }
      if (vm.count("campaigns")) {
          config->campaigns = vm["campaigns"].as<unsigned long>();
      }
      if (vm.count("cores")) {
          config->cores = vm["cores"].as<unsigned long>();
      } else { /* default value for # cores (save one for source) */
          if (config->cores <= 0 || config->cores >= std::thread::hardware_concurrency()) {
              config->cores = std::thread::hardware_concurrency() - 1;
          }
      }

      if (vm.count("input_file")) {
          config->input_file = vm["input_file"].as<vector<string>>()[0];
      }
  }
  catch(exception& e) {
      cerr << "error: " << e.what() << "\n";
      abort();
  }
  catch(...) {
      cerr << "Exception of unknown type!\n";
      abort();
  }
}
