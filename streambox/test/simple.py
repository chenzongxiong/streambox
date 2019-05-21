all_tests = [
	{
		"name" : "grep-10ms",  # this is very strict test
		"exec" : "./test-grep.bin", 
		"cores" : [56 - 1],
		"records" : 1000,
		"record_size" : 2000,
		"target_ms" : 10,
		"input_file" : "/ssd/1g.txt",
		# --- optional --- #
		"tput_baseline" : 375, # used to be compared with the test results		
		"tput_hint" : 400, 
		# --- control --- #
		"disable" : True	# skip the test		
	},
	{
		"name" : "grep-500ms",  
		"exec" : "./test-grep.bin", 
		"cores" : [56 - 1],
		"records" : 1000,
		"record_size" : 2000,
		"target_ms" : 500,
		"input_file" : "/ssd/1g.txt",
		# --- optional --- #
		"tput_baseline" : 375, # used to be compared with the test results		
		"tput_hint" : 4000, 
		# --- control --- #
		"disable" : True	# skip the test		
	},	
	{
		"name" : "grep-normal-500ms",  
		"exec" : "./test-grep.bin", 
		"cores" : 54,
		"records" : 1000 * 1000,
		"record_size" : 2000,
		"target_ms" : 500,
		"input_file" : "/ssd/9g.txt",
		# --- optional --- #
		"tput_baseline" : 10742, # used to be compared with the test results		
		"tput_hint" : 10000, 
		# --- control --- #
		#"disable" : True	# skip the test		
	},		
	{
		"name" : "wc-1sec",
		"exec" : "./test-wc.bin", 
		#"cores" : [56 - 1],
		"records" : 1000 * 1000,		# records per epoch
		"record_size" : 100,	
		"target_ms" : 1000,		
		"input_file" : "/ssd/1g.txt",
		# --- optional --- #
		"tput_baseline" : 4000, # used to be compared with the test results
		"tput_hint" : 4000,  # the throughput value that test should try first
		# --- control --- #
		"disable" : True	# skip the test
	},
]	

