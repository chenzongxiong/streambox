#include "papi_sampler.hpp"
using namespace std;

void handle_error(int retVal, std::string msg)
{
	std::cout << msg << std::endl;
	printf("PAPI error %d: %s\n", retVal, PAPI_strerror(retVal));
	exit(1);
}
void handle_error (int retval)
{
     printf("PAPI error %d: %s\n", retval, PAPI_strerror(retval));
     exit(1);
}

void PapiSampler::resetSampling()
{
//	int retVal = PAPI_cleanup_eventset(EventCodes);
//	if (retVal!= PAPI_OK)
//		handle_error(retVal, "Error for resetting counter");
//
//	retVal = PAPI_reset(EventCodes);
//	if (retVal!= PAPI_OK)
//		handle_error(retVal, "Error for resetting counter");
}

PapiSampler::PapiSampler(std::string& pSeq, std::string& pConfPath)
{
	confPath = pConfPath;
	counterSeq = pSeq;
	cout << "using papiFile=" << confPath << " preset=" << counterSeq << endl;
}

void PapiSampler::printSampling()
{
	std::stringstream ss_ret;
	long_long* overallRes = new long_long[numberOfCounters]{0};

	//create result string
	ss_ret << "Counter Names: ";
	for(size_t i = 0; i < numberOfCounters; i++)
	{
		ss_ret << counter[i] << ",";
	}
	ss_ret << std::endl;

	ss_ret << "Counter Values: " << endl;
	for(size_t th = 0; th < threadCnt; th++)
	{
		ss_ret << "\t Thread " << th << " ";
		for(size_t i = 0; i < numberOfCounters; i++)
		{
			ss_ret << result[th][i] << ",";
			overallRes[i] += result[th][i];
		}
		ss_ret << endl;
	}
	ss_ret << "CounterOverall=";
	for(size_t i = 0; i < numberOfCounters; i++)
	{
		ss_ret << overallRes[i]  << ",";
	}
	ss_ret << std::endl;
	std::cout << ss_ret.str();
}

//long_long* getSamplingResult()
//{
//
//}


bool PapiSampler::init(size_t pThreadCnt)
{
	std::ifstream infile(confPath);
	assert(infile.is_open());
	threadCnt = pThreadCnt;

	//find line in file
	std::string line = "";
	std::string delimiter = ",";
	std::string values;

	while (std::getline(infile, line))
	{
		 std::istringstream iss(line);
		 if(line.substr(0,counterSeq.size())== counterSeq)
		 {
			 size_t pos = line.find('=');
			 values = line.substr(pos+1, line.length());
			 break;
		 }

	}

	//split in counter values
	size_t pos = 0;
	std::string token;
	while ((pos = values.find(delimiter)) != std::string::npos) {
		token = values.substr(0, pos);
//	    std::cout << token << std::endl;
		counter.push_back(token);
		values.erase(0, pos + delimiter.length());
	}
	counter.push_back(values);

	//init events
	EventCodes = new int[counter.size()];
	result = new long_long*[threadCnt];
	for(size_t i = 0; i < threadCnt; i++)
	{
		result[i] = new long_long[counter.size()];
	}
	numberOfCounters = counter.size();

	//	/* Initialize the library */
	int retval = PAPI_library_init(PAPI_VER_CURRENT);
	if (retval != PAPI_VER_CURRENT) {
		printf("PAPI library init error!\n");
		exit(1);
	}

	if (PAPI_thread_init(pthread_self) != PAPI_OK)
	{
		printf("PAPI Thread init error!\n");
		 handle_error(1);
	}

	//fill array
	int i = 0;
	int eventCode;
	PAPI_event_info_t info;
	stringstream cntcout;
#ifndef SILENCE
	cntcout << counter.size() << " Counters: ";
#endif
	for (std::vector<std::string>::iterator it = counter.begin() ; it != counter.end(); ++it)
	{
		char* tempCode = const_cast<char*>((*it).c_str());
		retval =  PAPI_event_name_to_code( tempCode , &eventCode);
		if(retval != PAPI_OK)
			handle_error(retval, "error converting event ");

		int retVal = PAPI_get_event_info(eventCode, &info);
		if (retVal != PAPI_OK)
			 handle_error(retVal, "Error for event info");

		//Check if event is available
		retVal = PAPI_query_event(eventCode);
		if (retVal != PAPI_OK)
			handle_error(retVal, "No counter: " +  std::string(info.symbol));

		EventCodes[i] = eventCode;
		cntcout << info.symbol << ",";
		i++;
	}
	std::cout << cntcout.str() << std::endl;
	return true;
}

bool PapiSampler::stopSampling(size_t threadNo)
{
	//std::cout << "stopping: " << getId(id) << std::endl;
	int retVal = PAPI_stop_counters(result[threadNo], numberOfCounters);

	if (retVal != PAPI_OK)
		handle_error(retVal, "Error for stopping counter");

	return true;
}

bool PapiSampler::startSampling()
{
//	int* Event = eventSets[getId(id)];
	int retVal = PAPI_start_counters(EventCodes, numberOfCounters);

	if (retVal!= PAPI_OK)
		handle_error(retVal, "Error for starting counter");

	return true;
}
