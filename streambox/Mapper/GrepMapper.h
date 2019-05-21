
/* NB: This works on *unwindowed* records.
 *
 * output: the string that matches the regex
 *
 * for regex, see http://www.cplusplus.com/reference/regex/match_results/str/
 *
 * since this transform saves regex, it becomes stateful. can we make do_map()
 * static? however, this seems less compelling as each do_map() is doing
 * more substantial work than a virtual function lookup */

#ifndef GREP_MAPPER_H
#define GREP_MAPPER_H

#include <regex>
#include <re2/re2.h>

#include "Mapper/Mapper.h"
#include "Values.h"
#include <pcre.h>
using namespace std;

template <class InputT = string_range,
		class OutputT = string,
		template<class> class BundleT_ = RecordBundle>
class GrepMapper : public Mapper<InputT> {
//  using OutputBundleT = RecordBitmapBundle<OutputT>;
	using OutputBundleT = BundleT_<OutputT>;

public:
  GrepMapper(const char * str, string sname = "grep_mapper")
  : Mapper<InputT>(sname), ex(str), re(str), regex_str(str) {
  	assert(re.ok());
	
	//hym: for PCRE
	//const char *aStrRegex = regex_str.c_str();
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


  /* using pcre
   * return: # of emitted records
   */
  uint64_t do_map_pcre(Record<InputT> const & in,
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
				output_bundle->add_record(Record<OutputT>(
								//psubStrMatchStr,
								//string(p, regex_str.size()
								string(testStrings + subStrVec[2*i], subStrVec[2*i+1] - subStrVec[2*i]),
								in.ts)
				);
				start_offset = subStrVec[2*i+1]; //XXX moved in ??
			}	
			//start_offset = subStrVec[1]; //should be moved into for?? XXX
		}
	}
out:
	return count;

  }


  /* using std::regex
   * return: # of emitted records */
//  uint64_t do_map_std(Record<InputT> const & in,
  uint64_t do_map_std(Record<InputT> const & in,
  		shared_ptr<OutputBundleT> output_bundle) {
	  	using namespace std::regex_constants;
	  	uint64_t count = 0;

			/* dirty hacking: we force 0-terminating each record's
 			 * string range. this is because that std::regex
			 * only works with 0-terminated string.
			 */
			*((char *)in.data.data + in.data.len - 1) = '\0';

			cregex_iterator it(in.data.data,
					in.data.data + in.data.len - 1,  /* be very careful here */
					ex);
			cregex_iterator reg_end;

			for (; it != reg_end; ++it) {
  			output_bundle->add_record(Record<OutputT>(
  					it->str(), /* XXX optimize to avoid copy? */
  					in.ts)
  			);
  			count ++;
//  			cout << it->str();
  		}
  	return count;
  }

  /* the version using Google's RE2
   * ref: https://gist.github.com/chezou/1395527
   *
   * return: # of emitted records */
  uint64_t do_map_re2(Record<InputT> const & in,
  		shared_ptr<OutputBundleT> output_bundle) {

			/* dirty hacking: we force 0-terminating each record's
 			 * string range. this is because that std::regex
			 * only works with 0-terminated string.
			 */
			*((char *)in.data.data + in.data.len - 1) = '\0';

			string matched;

#if 0	/* only match once in a record */
			if (RE2::PartialMatch(in.data.data, re, &matched)) {
				/* XXX how do we deal with multiple matches in one record?
				 * use StringPiece? */
				output_bundle->add_record(Record<OutputT>(matched, in.ts));
		  	return 1;
			} else
				return 0;
#endif

			/* using a local re -- even slower */
//			RE2 re1(regex_str);
//			assert(re1.ok());

#if 1
			re2::StringPiece input(in.data.data);
			uint64_t count = 0;
			while(RE2::FindAndConsume(&input, re, &matched)) {
				output_bundle->add_record(Record<OutputT>(matched, in.ts));
				count ++;
			}
			return count;
#endif
  }

  uint64_t do_map_strstr(Record<InputT> const & in,
  		shared_ptr<OutputBundleT> output_bundle) {

		*((char *)in.data.data + in.data.len - 1) = '\0';

		uint64_t count = 0;
//		char *p = strstr((char *)in.data.data, regex_str.c_str());
		char *p = (char *)in.data.data;

		while ((p = strstr(p, regex_str.c_str()))) {
			count ++;
			output_bundle->add_record(Record<OutputT>(string(p, regex_str.size()), in.ts));
			p += regex_str.size();
			assert(p <= (char *)in.data.data + in.data.len - 1);
		}

		return count;
  }


  uint64_t do_map(Record<InputT> const & in,
  		shared_ptr<OutputBundleT> output_bundle) {
//  	return do_map_re2(in, output_bundle);
//  	return do_map_std(in, output_bundle);
//  	return do_map_strstr(in, output_bundle);
  	return do_map_pcre(in, output_bundle);
  }

  void ExecEvaluator(int nodeid, EvaluationBundleContext *c,
  		shared_ptr<BundleBase> bundle_ptr = nullptr) override;

private:
  std::regex ex;
  RE2 re;

  //hym: for PCRE
  pcre *reCompiled;
  const char *pcreErrorStr;
  int pcreErrorOffset;
  pcre_extra *pcreExtra;

public:
  string regex_str;
};

#endif // GREP_MAPPER_H
