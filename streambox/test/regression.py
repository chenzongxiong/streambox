all_tests = [
	{
		"name" : "wc-fast",
		"exec" : "./test-wc.bin", 
		"records" : 1000 * 1000,		# records per epoch
		"record_size" : 100,	
		"target_ms" : 1000,		
		"input_file" : "/ssd/1g.txt",
		# --- optional: soft delay --- #
		#"softdelay_maxbad_ratio" : 0.1, # okay if aonmaly delay % is less than this in a window
		#"softdelay_maxbad_ms" : 2000, # upper bound of anonmaly delay
		# --- optional --- #
		# "cores" : 54,   # if unspecified, fall back to app default
		"tput_baseline" : 5200, # used to be compared with the test results
		"tput_hint" : 5000,  # the throughput value that test should try first
		# --- control --- #
#		"disable" : True	# skip the test
	},
	
	{
		"name" : "wingrep-fast",
		"exec" : "./test-wingrep.bin", 
		"records" : 1000 * 1000,		# records per epoch
		"record_size" : 1000,	
		"target_ms" : 1000,		
		"input_file" : "/ssd/9g.txt",
		# --- optional --- #
		# "cores" : 54,   # if unspecified, fall back to app default		
		"tput_baseline" : 38500, # used to be compared with the test results
		"tput_hint" : 37000,  # the throughput value that test should try first
		# --- control --- #
#		"disable" : True	# XXX skip the test
	},

	{
		"name" : "test-join-2-fast",
		"exec" : "./test-join-2.bin", 
		"records" : 1000 * 1000,		# records per epoch
		"record_size" : 8, #sizeof(long)	
		"target_ms" : 1000,		
		"input_file" : "/ssd/test-digit.txt",
		# --- optional --- #
		# "cores" : 54,   # if unspecified, fall back to app default		
		"tput_baseline" : 5200, # used to be compared with the test results
		"tput_hint" : 5000,  # the throughput value that test should try first
		# --- control --- #
#		"disable" : True	# XXX skip the test
	},
	
	{
		"name" : "test-distinct-fast",
		"exec" : "./test-distinct.bin", 
		"records" : 1000 * 1000,		# records per epoch
		"record_size" : 100,
		"target_ms" : 1000,		
		"input_file" : "/ssd/train-out.txt",
		# --- optional --- #
		# "cores" : 54,   # if unspecified, fall back to app default		
		"tput_baseline" : 2000, # xzl: can do 2000? used to be compared with the test results
		"tput_hint" : 2000,  # the throughput value that test should try first
		# --- control --- #
#		"disable" : True	# XXX skip the test
	},
	
	{
		"name" : "networklatency-fast",
		"exec" : "./networklatency.bin", 
		"records" : 500 * 1000,		# records per epoch
		"record_size" : 40, #sizeof(struct srcdst_rtt)
		"target_ms" : 1000,		
		"input_file" : "/ssd/region_Raw_PingmeshData.result",
		# --- optional --- #
		# "cores" : 54,   # if unspecified, fall back to app default		
		"tput_baseline" : 1000, # xzl: 878
		"tput_hint" : 800,  # the throughput value that test should try first
		# --- control --- #
#		"disable" : True	# XXX skip the test
	},
	
	{
		"name" : "test-tweet-fast",
		"exec" : "./test-tweet.bin", 
		"records" : 1000 * 1000,		# records per epoch
		"record_size" : 200,
		"target_ms" : 1000,		
		"input_file" : "/ssd/twitter_download/filtered_tweets.txt",
		# --- optional --- #
		# "cores" : 54,   # if unspecified, fall back to app default		
		"tput_baseline" : 5000, # used to be compared with the test results
		"tput_hint" : 4000,  # the throughput value that test should try first
		# --- control --- #
#		"disable" : True	# XXX skip the test
	},
]	

