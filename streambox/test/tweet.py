all_tests = [	
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

