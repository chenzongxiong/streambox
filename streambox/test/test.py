#!/usr/bin/env python

# Copyright (c) 2011 The Chromium Authors. All rights reserved.
# Use of this source code is governed by a BSD-style license that can be
# found in the LICENSE file.

"""Android system-wide tracing utility.

This is a tool for capturing a trace that includes data from both userland and
the kernel.	It creates an HTML file for visualizing the trace.
"""

import errno, optparse, os, select, subprocess, sys, time, zlib, signal
import re
import tempfile
import math
import json

from regression import *
#from simple import *
#from simple1 import *		# for testing this tester per se

"""
test configuration
"""
config_tput_min = 100
config_tput_max = 50000
config_tput_resolution = 100
config_output_timeout_sec = 10  # can be floating pt
config_max_runtime_sec = 60 	# give up 

"""
default args
"""
#config_default_cores = [4, 12, 32, 56]
config_default_records = 1000 * 1000	# records per epoch

'''
global vars
'''
the_output_dir = ""

"""
test cases
"""

"""
	{
		"name" : "grep",
		"exec" : "./test-grep.bin", 
		#"cores" : [4, 12, 32, 56],
		"cores" : [56],
		"records" : 1000,
		"record_size" : 2000,
		"target_ms" : 1000,
		"input_file" : "/ssd/1g.txt",
		# --- optional --- #
		"tput_hint" : 4000, 
	},
"""

"""
app_list = [
		"test-grep", 
		"test-wc",
		"test-wingrep",
		"test-join",
		"test-join-2",
		"networklatency",
		"test-distinct",
		"test-tweet"
		]
"""

"""
sample line:
dump markers: >>>>>>>>>>>>>>>>>>>>>>>>>>>>total 7 ms
# return: delay in ms
"""
def parse_line_delay(line):
	delay_regex = r'''dump markers: >>>>>>>>>>>>>>>>>>>>>>>>>>>>total (\d+) ms'''
	m = re.match(delay_regex, line)
	if not m:
		return None
	else:
		return (m.group(1))

"""
sample line (for backward compatibility; future sources should have same name):
          unbounded-inmem     19.07     19.07     19.53     19.53
              [unbounded]     20.35     20.35   2604.17   2604.17
              [netsource]     31.79     31.79    813.80    813.80
              
return: tput in krec/s (floating)
"""          
def parse_line_tput(line):
	# the "?:" is for "non capturing" group. 
	regex = r'''\s+(?:unbounded-inmem|\[unbounded\]|\[netsource\])\s+(\d+\.\d+)\s+(\d+\.\d+)\s+(\d+\.\d+)\s+(\d+\.\d+)'''
	m = re.match(regex, line)
	if not m:
		return None
	else:
		recent_mbs = m.group(1)
		avg_mbs = m.group(2)
		recent_krecs = m.group(3)
		avg_krecs = m.group(4)
		return (avg_krecs)

# stateless 
# @delays: a list of all historical delays
# return: 

DECIDE_FAIL 	= 1		# failed. abort.
DECIDE_CONT 	= 2		# should continue
DECIDE_OK 	= 3		# target_delay is met
#DECIDE_DUNNO 	= 4		# can't decide yet
DECIDE_EXCEPTION 	= 5		# what happened?

decide_descs = ["", "fail", "cont", "ok", "dunno", "exception"]

# XXX catch c-c signal to ensure all test programs killed XXX

def decide_delay(delays, target_delay, test):
	n_init_samples = 10 # drop them all. must >= 1 since we do i-1 below  
	n_recent_samples = 15 
	n_total_samples = len(delays)
	
	# must at least pass the init phase and have 1 window
	if len(delays) < n_init_samples + n_recent_samples: 
		return DECIDE_CONT
		
	# check most recent N samples 
	is_go_up = True
	is_go_down = True
	n_good = 0
	n_bad = 0
	n_anomaly = 0
	
	# --- decide by trend --- #
	# trend -- do we keep going up/down?  
	# compute mov avg windows
	n_moving_avg_win = 3
	assert(n_recent_samples - n_moving_avg_win >= 5)
	
	if (n_recent_samples < n_moving_avg_win): # consider increasing n_recent_samples
		return 	DECIDE_CONT
	
	mov_delays = []	
	for i in range(n_total_samples - n_recent_samples, n_total_samples - n_moving_avg_win):
		s = 0
		for offset in range(n_moving_avg_win):
			s += delays[i + offset]
		s /= n_moving_avg_win
		mov_delays.append(s)

	if (len(mov_delays) <= 1):
		return 	DECIDE_CONT
		 
	for i in range(len(mov_delays) - 1):
		if mov_delays[i] >= mov_delays[i+1]: # since we're mv avg, robust to ==
			is_go_up = False
		if mov_delays[i] <= mov_delays[i+1]: # since we're mv avg, robust to ==
			is_go_down = False
			
	#if (is_go_up and n_good > n_recent_samples / 2):
	if (is_go_up):  	# shall we also say how far it is from our target delay?
		print "latency fail: rising delays", delays[n_total_samples - n_recent_samples:]
		return DECIDE_FAIL

	# --- decide by target delay --- #
	# all good (or good enough for softdelay) -- okay
	# all bad -- fail 
	# a mix of good & bad -- undecided
	#elif (1.0 * n_good / n_recent_samples <= 0.05)
	maxgood = target_delay
	maxbad = target_delay # softdelay can't pass if any sample larger than this. but may not fail immediately
	max_n_bad = 0  # pass , if #bads in a window smaller than this
	if test.has_key("softdelay_maxbad_ratio"):
		max_n_bad = math.ceil(n_recent_samples * test["softdelay_maxbad_ratio"])
		maxbad = test["softdelay_maxbad_ms"]

	max_sofar = -1 

	# go through individual samples
	for i in range(n_total_samples - n_recent_samples, n_total_samples):
		if delays[i] > maxgood:
			#print "over target", i, delays[i]
			n_bad += 1
		if delays[i] > max_sofar:  
			#print "latency fail. anonmaly:", delays[i], "over limit", maxbad
			#return DECIDE_FAIL
			max_sofar = delays[i] 

	if (n_bad <= max_n_bad and max_sofar < maxbad): # good (enough)
		print "latency okay! bad:", n_bad, delays[n_total_samples - n_recent_samples:]
		return DECIDE_OK

	if (n_bad == n_recent_samples):  # all bad
		print "latency fail: all over targets", delays[n_total_samples - n_recent_samples:]
		return DECIDE_FAIL

	# undecided -- a mix of good and bad
	print "can't decide. n_bad", n_bad, "maxsofar", max_sofar, "maxbad", maxbad, \
		   "recent delays: ", delays[n_total_samples - n_recent_samples:]

	return 	DECIDE_CONT

	# XXX more logic? XXX
	
	#elif (is_go_down)
	#	return DECIDE_CONT

is_stop = False
def stop_test_handler(signal, frame):
	print "Stopping the test.... (ctrl-c)"
	is_stop = True

'''
return (status, tput)
tput is floating pt. <0 if none achieved
'''
def run_one_test(test, atrace_args):
	trace_started = False
	leftovers = ''
	start_sec = time.time()
	
	#html_filename = os.path.join(the_output_dir, "%s-tput%d.log" %(test["name"], test["target_tput"]))
	# will rename this after test done
	html_filename = os.path.join(the_output_dir, test["name"], 
					"%s-tput%d-ongoing.log" %(test["name"], test["target_tput"]))
	html_file = open(html_filename, 'w')
	
	delays = []
	avg_tput = -1
	
	test_status = DECIDE_EXCEPTION
	
	# xzl -- run the actual command on target --
	print "run cmd ===>", atrace_args
	adb = subprocess.Popen(atrace_args, stdout=subprocess.PIPE,
								 stderr=subprocess.PIPE, 
								 stdin=subprocess.PIPE)   # to feed stdin

	# per Python doc, this may deadlock; but it is NON BLOCKING!
	adb.stdin.write('asdb\r\nsdfsd\rdsfsdf\nsfds\r\ndsfsd\r\n')
	# this will block read stdout/stderr... bad
	#adb.communicate(input=b'asdb\r\nsdfsd\rdsfsdf\nsfds\r\ndsfsd\r\n')[0]
	
	# XXX catch expcetion to kill process
	while True:
		# we may interrupt the blocking call. 
		# http://stackoverflow.com/questions/5633067/signal-handling-in-pylons
		try:
			ready = select.select([adb.stdout, adb.stderr], [], [adb.stdout, adb.stderr],
							 config_output_timeout_sec)
		except select.error, v:
			if v[0] != errno.EINTR: 
				print "don't know how to handle"
				raise
			else: 
				print "intr!"
		
		if adb.stderr in ready[0]:
			err = os.read(adb.stderr.fileno(), 4096)
			sys.stderr.write(err)
			sys.stderr.flush()
		if adb.stdout in ready[0]:
			#print "got one std line"
			out = leftovers + os.read(adb.stdout.fileno(), 4096)
			#if options.from_file is None:
			if True: # xzl
				out = out.replace('\r\n', '\n')
			if out.endswith('\r'):
				out = out[:-1]
				leftovers = '\r'
			else:
				leftovers = ''

			# XXX toggle this 
			#sys.stdout.write(out)
			#sys.stdout.flush()
			
			html_file.write(out)
			html_file.flush()
			
			lines = out.splitlines(True)
			out = ''
			for i, line in enumerate(lines):
			
				tput = parse_line_tput(line)
				if tput != None:
					avg_tput = tput
					#print "XXXXXXXXXX Got tput", avg_tput
					#sys.exit(-1)
					
				# ---- xzl: parse an output line ---- #
				#sys.stdout.write("collecting trace...")
				#sys.stdout.flush()
				delay = parse_line_delay(line)
				
				
				if delay != None:
					delays += [int(delay)]
					decision = decide_delay(delays, int(test["target_ms"]), test)
					if (decision == DECIDE_OK and avg_tput > 0):
						# target delay met, we're done
						print "test-wrapper: okay to meet latency; stop"
						adb.kill()
						test_status = DECIDE_OK
						break;  # will go check status
					elif (decision == DECIDE_FAIL and avg_tput > 0):
						print "test-wrapper: fail to meet latency; stop"
						adb.kill()
						test_status = DECIDE_FAIL
						break;  # will go check status
				#out = ''.join(lines[i:]) #xzl: don't skip any line
				#html_file = open(html_filename, 'w')
				#html_file.write("# " + "%s" %atrace_args)	#xzl: save our command for dbg
				trace_started = True
				
			#sys.stdout.write(out)
			#sys.stdout.flush()
			
			#html_out = out.replace('\n', '\\n\\\n')
			#html_out = out
			#if len(html_out) > 0:
			#	html_file.write(html_out)
				
		# xzl -- done reading a wave of output from target -- 

		if (time.time() - start_sec > config_max_runtime_sec):
			print >> sys.stderr, "test timeout. after %.2f sec still can't decide tput" %(time.time() - start_sec)
			adb.kill()
			test_status = DECIDE_FAIL  # it's not an expcetion: maybe the latency just can't stablize.
			time.sleep(1)

		# check prog status
		result = adb.poll()
		if result is not None:
			break
		# result == None means child not end yet
		
	if result != 0:
		#print >> sys.stderr, 'program returned error code %d' % result
		#print >> sys.stderr, result
		print >> sys.stderr, 'program killed'
		pass
	elif trace_started:	 # xzl: program ends okay and we've collected some trace
		html_out = dec.flush().replace('\n', '\\n\\\n').replace('\r', '')
		if len(html_out) > 0:
			html_file.write(html_out)
		#html_file.write(html_suffix)
		html_file.close()
		#print " done\n\n		wrote file://%s/%s\n" % (os.getcwd(), options.output_file)
		print " done\n\n		wrote %s\n" % (options.output_file)		
	else:
		print >> sys.stderr, ('An error occured while capturing the trace.	Output ' +
			'file was not written.')
	
	html_file.close()
	
	# rename based on results
	
	if (test_status == DECIDE_FAIL):
		nname = os.path.join(the_output_dir, test["name"], 
				"%s-fail-target_tput%d.log" %(test["name"], test["target_tput"]))
	elif (test_status == DECIDE_EXCEPTION):
		nname = os.path.join(the_output_dir, test["name"], 
				"%s-exception-target_tput%d.log" %(test["name"], test["target_tput"]))
	elif (test_status == DECIDE_OK):
		nname = os.path.join(the_output_dir, test["name"], 
				"%s-ok-target_tput%d-actual%d.log" %(test["name"], test["target_tput"], int(float(avg_tput))))
	else:
		assert(False)

	os.rename(html_filename, nname)
	
	return test_status, float(avg_tput)
	
def save_test_config(test):
	config_fname = os.path.join(the_output_dir, test["name"], "config.json")
	cf = file(config_fname, "w")
	# pretty, see https://docs.python.org/2/library/json.html
	json.dump(test, cf, indent=4, separators=(',', ': ')) 
	cf.close()

# @core == -1 if unspecified on cmdline
def print_test_cmd(test, core, tput_best):
	print "---------------------"
	print "you can repeat the best experiment by:"
	print os.path.abspath(test["exec"]), "\\ \n", "--target_tput", tput_best, \
			"--records", test["records"], "\\ \n", \
			"--input_file", test["input_file"], \
			"--record_size", test["record_size"],   # add more?
	if (core != -1):
		print "--cores", core
	else:
		print 
	print "---------------------"

# return: actual_tput
# XXX only support one core now. but that's fine
def launch_one_test(test):
	tput = 1000 # the target_tput passed to cm
	tput_best = -1 # the target_tput corresponding to the best actual tput
	tput_low = config_tput_min 
	tput_high = config_tput_max
	
	if test.has_key("disable") and test["disable"]:
		print >> sys.stderr, "skip disabled test: ", test["name"]
		return
		
	if test.has_key("tput_hint"):
		tput = test["tput_hint"]
	
	os.mkdir(os.path.join(the_output_dir, test["name"]))
	
	save_test_config(test)

	start_sec = time.time()
	
	core = -1
	if test.has_key("cores"):
		core = test["cores"]

	actual_tput = -1	
	while True: # execute the test with different tput
	
		# xzl: the actual command. all args except tput	
		print "---------> run with tput %d (low %d high %d) target %d ms" %(tput, tput_low, tput_high, test["target_ms"])
		
		args = [test["exec"], 
					#'--cores=%s' %core,
					'--target_tput=%d' %tput,
					'--records=%s' %test["records"],
					'--input_file=%s' %test["input_file"],
					'--record_size=%s' %test["record_size"],
					# todo: add more #
					]
		if (core != -1):
			args.append('--cores=%d' %core)
			
		test["target_tput"] = tput  # the runner can gen filename based on this
		
		status, t = run_one_test(test, args)			
		
		if (status == DECIDE_EXCEPTION):
			# save test exec time
			test["elapsed_sec"] = time.time() - start_sec
			return -1
		if (status == DECIDE_OK):			
			if (t > 0): # update actual tput
				print "actual_tput prev:", actual_tput, "new:", t
				#assert(actual_tput < 0 or t >= actual_tput + config_tput_resolution) # this may happen...?
				#assert(actual_tput < 0 or 1.05 * t >= actual_tput) # this may happen...?			
				if not (actual_tput < 0 or 1.05 * t >= actual_tput): 
					# after lowering target tput, we achieve target lat, 
					# but the actual tput is actually lower than prev. we're done.
					test["elapsed_sec"] = time.time() - start_sec
					print_test_cmd(test, core, tput_best)
					return actual_tput
				else:
					actual_tput = t
					tput_best = tput
			# --- is range small enough? --- #
			if (tput_low + config_tput_resolution > tput_high): # range small enough, done.
				print_test_cmd(test, core, tput_best)
				# save test exec time
				test["elapsed_sec"] = time.time() - start_sec
				return actual_tput
			else:
				# we set tput based on our prior tput setting, not the actual @actual_tput
				tput_low = tput
				tput = (tput_low + tput_high) / 2 
				continue
		elif (status == DECIDE_FAIL):
			if (tput_low + config_tput_resolution > tput_high):
				# save test exec time
				test["elapsed_sec"] = time.time() - start_sec				
				return actual_tput  # == -1 if we never seen one
			else:
				tput_high = tput
				tput = (tput_low + tput_high) / 2
				continue				
		else:
			print "err?"
			sys.exit(1)
	
		
if __name__ == '__main__':
	signal.signal(signal.SIGINT, stop_test_handler)
	
	the_output_dir = tempfile.mkdtemp()
	
	# results will be inserted in place into @all_tests
	
	''' check & print all test info '''
	test_names = {}
	# detect duplicate test names 
	for test in all_tests:
		if test_names.has_key(test["name"]):
			print >> sys.stderr, "abort: duplicate test names:", test["name"];
			sys.exit(1)
		test_names[test["name"]] = 1  # remember to check duplicates
	
		if test.has_key("softdelay_maxbad_ms") and test["softdelay_maxbad_ms"] < test["target_ms"]:
			print >> sys.stderr, "abort: config err: [%s] softdelay maxbad ms < target ms" %test["name"]
			sys.exit(-1)

	''' print menu '''
	print "========================================"
	print "select tests to run (enter to run all)"
	for i, test in enumerate(all_tests):
		print i, test["name"];
	try:
		choice = int(raw_input('Enter your input:'))
		print "Okay ... Will run test", all_tests[choice]["name"]
		all_tests = [all_tests[choice]]
	except ValueError:
		print "Okay ... Will run all tests."

	for test in all_tests:
		atput = launch_one_test(test)
		if atput < 0:
			print "%s exception: can't get the tput." %test["name"]
			test["actual_tput"] = -1  # is this okay?
		else:
			test["actual_tput"] = atput
			print "%s completes: actual tput %d krec/s target_delay %d ms" %(test["name"], atput, test["target_ms"])
	
	print "========================================"
	print "%20s %10s %10s %10s %6s %15s" %("test", "target_ms", "tput/krecs", "base/krecs", "improve%", "elapsed/sec")	
	for test in all_tests:
		tput_inc = -999.99
		tput_inc_str = "--"
		tput_baseline_str = "--"
		
		if test.has_key("disable") and test["disable"]:
		#if not test.has_key("elapsed_sec"): # test never executed?
			print "%10s -- skipped -- " %(test["name"])
			continue
		if test.has_key("tput_baseline"):
			tput_inc = 100.0 * (test["actual_tput"] - test["tput_baseline"]) / test["tput_baseline"]
			tput_inc_str = "%.2f" %(tput_inc)
			tput_baseline_str = "%d" %(test["tput_baseline"])
			#print "baseline is", test["tput_baseline"]
			
		print "%20s %10d %10d %10s %6s %15.2f"  \
				%(test["name"], test["target_ms"], test["actual_tput"], tput_baseline_str, tput_inc_str, test["elapsed_sec"])
	print "========================================"
	print "diff=-999 means no baseline provided"
	
	print "all done. check result dir:\n ls ", the_output_dir


