#ifndef WIN_GREP_MAPPER_H
#define WIN_GREP_MAPPER_H

#if 0
#define DECL_EVAL2(T) \
	using TransformT = T<InputT,OutputT>; \
	using EvalT = T##Evaluator<InputT,OutputT>

template<typename I,typename O> class WindowedGrepMapperEvaluator;
class EvaluationBundleContext;
#endif

// this works on windowed records.
//
// stateless
#include <string>
#include <regex>
#include <re2/re2.h>

#include "Mapper/Mapper.h"
#include "Values.h"
#include <pcre.h>
using namespace std;

template <class InputT = string_range,
		typename OutputT = creek::string>
class WindowedGrepMapper : public Mapper<InputT> {
/*	DECL_EVAL2(WindowedGrepMapper); */
	using OutputBundleT = WindowsBundle<OutputT>;

public:
  WindowedGrepMapper(const char * str, string sname = "grep_mapper")
  : Mapper<InputT>(sname), ex(str), re(str), regex_str(str) { 
  	
	  //hym: for PCRE
	  reCompiled = pcre_compile(str, 0, &pcreErrorStr, &pcreErrorOffset, NULL);
	  if(reCompiled == NULL) {
		  std::cout << "ERROR: Could not compile PCRE: " << regex_str << std::endl;
		  exit(1);	 
	  }	

	  pcreExtra = pcre_study(reCompiled, 0, &pcreErrorStr);
	  if(pcreErrorStr != NULL) {
		  std::cout << "ERROR: Could not study " << str << " " << pcreErrorStr << std::endl;
		  exit(1);
	  }
  }

 // return: # of emitted records 
  uint64_t do_map_pcre(Window const & win, Record<InputT> const & in,
  		shared_ptr<OutputBundleT> output_bundle) {

	  uint64_t count = 0;

	  /* dirty hacking: we force 0-terminating each record's
	   * string range. this is because that std::regex
	   * only works with 0-terminated string.
	   */
	  *((char *)in.data.data + in.data.len - 1) = '\0';

	  int pcreExecRet;
	  int subStrVec[30];
	  int i;    
	  char * testStrings = (char *)in.data.data;
	  int subject_length = in.data.len;
	  int start_offset = 0;
	  while(1){
		  pcreExecRet = pcre_exec(reCompiled,
				  pcreExtra,
				  testStrings,
				  subject_length,
				  start_offset, //start looking at this point
				  0, 	      //OPTIONS
				  subStrVec,    //subStrVec[0] is the position of the first matched substring
				  30);
		  if(pcreExecRet < 0) {
			  //no match
			  goto out;
		  }else{
			  count += pcreExecRet;
			  for (i = 0; i < pcreExecRet; i++){
				  //subject_length = subject_length - (subStrVec[2*i+1] - subStrVec[2*i]);
				  //pcre_get_substring(testStrings, subStrVec, pcreExecRet, i, &(psubStrMatchStr));
				  //std::cout << "matched string: " << psubStrMatchStr << std::endl;
				  output_bundle->add_value(win, Record<OutputT>(
							  //psubStrMatchStr, /* XXX optimize to avoid copy? */
							  creek::string(testStrings + subStrVec[2*i], subStrVec[2*i+1] - subStrVec[2*i]),
							  in.ts));
				  start_offset = subStrVec[2*i+1]; //XXX moved in ??

			  }	
			  //start_offset = subStrVec[1]; //XXX should be moved into for loop??
		  }
	  }
out:
	  return count;
  }

  // return: # of emitted records
  uint64_t do_map_std(Window const & win, Record<InputT> const & in,
  		shared_ptr<OutputBundleT> output_bundle) {
	  	using namespace std::regex_constants;
	  	uint64_t count = 0;

		 // dirty hacking: we force 0-terminating each record's
			// string range. this is because that std::regex
			// only works with 0-terminated string.

		 *((char *)in.data.data + in.data.len - 1) = '\0';

			cregex_iterator it(in.data.data,
					in.data.data + in.data.len - 1,  /* be very careful here */
					ex);
			cregex_iterator reg_end;

			for (; it != reg_end; ++it) {
  			output_bundle->add_value(win, Record<OutputT>(
  					it->str(), /* XXX optimize to avoid copy? */
  					in.ts));
  			count ++;
  		}

  	return count;
  }

  /* re2 */
  uint64_t do_map_re2(Window const & win, Record<InputT> const & in,
  		shared_ptr<OutputBundleT> output_bundle) {

			/* dirty hacking: we force 0-terminating each record's
 			 * string range. this is because that std::regex
			 * only works with 0-terminated string.
			 */
			*((char *)in.data.data + in.data.len - 1) = '\0';

			string matched;

			re2::StringPiece input(in.data.data);
			uint64_t count = 0;
			while(RE2::FindAndConsume(&input, re, &matched)) {
				output_bundle->add_value(win, Record<OutputT>(creek::string(matched), in.ts));
				count ++;
			}
			return count;
  }

  uint64_t do_map_strstr(Window const & win, Record<InputT> const & in,
  		shared_ptr<OutputBundleT> output_bundle) {

		*((char *)in.data.data + in.data.len - 1) = '\0';

		uint64_t count = 0;
		char *p = (char *)in.data.data;

		while ((p = strstr(p, regex_str.c_str()))) {
			count ++;
			output_bundle->add_value(win, Record<OutputT>(string(p, regex_str.size()), in.ts));
			p += regex_str.size();
			assert(p <= (char *)in.data.data + in.data.len - 1);
		}

		return count;
  }

  uint64_t do_map(Window const & win, Record<InputT> const & in,
  		shared_ptr<OutputBundleT> output_bundle) {
//  	return do_map_std(win, in, output_bundle);
//  	return do_map_re2(win, in, output_bundle);
//  	return do_map_pcre(win, in, output_bundle);
			return do_map_strstr(win, in, output_bundle);
  }

  /* instantiate a corresponding evaulator and run it (in sync) */
#if 0
  void RunEvaluator(int nodeid, PTransform *transform,
  		EvaluationBundleContext *c) {
  			EvalT eval(nodeid);
  			eval.evaluate(dynamic_cast<TransformT*>(transform), c);
  }

  void ExecEvaluator(int nodeid, EvaluationBundleContext *c) override {
		EvalT eval(nodeid); /* instantiate an evaluator */
		eval.evaluate(this, c);
  }
#endif

  void ExecEvaluator(int nodeid, EvaluationBundleContext *c,
  		shared_ptr<BundleBase> ptr = nullptr) override;

private:
  std::regex ex;
  RE2 re;
 
  //hym: for PCRE
  pcre *reCompiled;
  pcre_extra *pcreExtra;
  const char *pcreErrorStr;
  int pcreErrorOffset;
public:
  string regex_str;
};

#endif // WIN_GREP_MAPPER_H
