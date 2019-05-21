all_tests = [
	{
		"name" : "wc-fast",
		"exec" : "./test-wc.bin", 
		#"cores" : [4, 12, 32, 56],  
		"cores" : 54,
		"records" : 1000,		# records per epoch
		"record_size" : 100,	
		"target_ms" : 1000,		
		"input_file" : "/ssd/1g.txt",
		# --- optional --- #
		"tput_baseline" : 5200, # used to be compared with the test results
		"tput_hint" : 5000,  # the throughput value that test should try first
		# --- control --- #
		#"disable" : True	# skip the test
	},
	{
		"name" : "wingrep-fast",
		"exec" : "./test-wingrep.bin", 
		"records" : 1000,		# records per epoch
		"record_size" : 1000,	
		"target_ms" : 1000,		
		"input_file" : "/ssd/9g.txt",
		# --- optional --- #
		"tput_baseline" : 38500, # used to be compared with the test results
		"tput_hint" : 37000,  # the throughput value that test should try first
		# --- control --- #
		#"disable" : True	# skip the test
	},
]
