#pragma once
#include <papi.h>
#include <fstream>
#include "assert.h"
#include <sstream>
#include <iostream>
#include <vector>

class PapiSampler
{
public:
	PapiSampler(std::string& seq, std::string& confPath);
	bool init(size_t threadCnt);
	bool startSampling();
	bool stopSampling(size_t threadNo);
	void printSampling();
	void resetSampling();
	long_long* getSamplingResult();

private:
	int* EventCodes;
	long_long** result;
	size_t numberOfCounters;
	std::string counterSeq;
	std::string confPath;
	std::vector<std::string> counter;
	size_t threadCnt;
};
