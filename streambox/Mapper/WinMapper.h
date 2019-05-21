#ifndef WINMAPPER_H_
#define WINMAPPER_H_

#include <fstream>

#include "Values.h"
#include "Mapper/Mapper.h"
#include "boost/date_time/posix_time/ptime.hpp"
#include <boost/algorithm/string.hpp>  /* split string */

template <class InputT, class OutputT, template<class> class BundleT>
class WinMapper : public Mapper<InputT> {

	using InputBundleT = BundleT<InputT>;
	using OutputBundleT = WindowsBundle<OutputT>;

public:
  const boost::posix_time::time_duration window_size;
  const ptime start; // the starting point of windowing.

private:
  const char *keyword_file_path = "/ssd/Twitter-Sentiment-Analysis/AFINN-111.txt";
  std::unordered_map<string, int> keywords;

public:
	WinMapper(string name,
      boost::posix_time::time_duration window_size,
      ptime start = Window::epoch)
	 : Mapper<InputT>(name), window_size(window_size), start(start) {

		/* load the keywords from file */
		int val, cnt = 0, mcnt = 0;
		std::ifstream infile(keyword_file_path);

#if 0 /* does not work since keywords may container space */
		while (infile >> word >> val)
		{
			cout << word << val << endl;
			cnt ++;
		}
#endif

		/* split strings. see
		 * http://stackoverflow.com/questions/236129/split-a-string-in-c
		 */
		std::vector<std::string> tokens;
		std::string line;
		while (std::getline(infile, line)) {
			boost::split(tokens, line, boost::is_any_of("\t "));

#if 0			/* debug */
			for(auto & s : tokens)
				cout << s << " ";
#endif

			int ntokens = tokens.size();
			xzl_bug_on(ntokens < 2);

			/* assemble the words from tokens... */
			string word = tokens[0];
			if (ntokens > 2) {
				for (int i = 1; i < ntokens - 1; i++) {
					word += string(" ") + tokens[i];
				}
				mcnt ++;
			}
//			cout << word << "*" << endl;
			val = stoi(tokens[ntokens - 1]);
//			cout << val << endl;
			cnt ++;

			keywords[word] = val;
		}

		if (cnt < 2000) {
			EE("there seems no enough keywords from %s", keyword_file_path);
			xzl_bug();
		} else
			EE("%d tweet keywords loaded from %s (%d has multi phrases and are ignored now)",
					cnt, keyword_file_path, mcnt);

	}


	void ExecEvaluator(int nodeid, EvaluationBundleContext *c,
	  		shared_ptr<BundleBase> bundle_ptr = nullptr) override;

//	uint64_t do_map(Window const & win, Record<InputT> const & in,
//	  		shared_ptr<OutputBundleT> output_bundle) {
//
//
//		output_bundle->add_value()
//	}

  uint64_t do_map(Record<InputT> const & in,
  		shared_ptr<OutputBundleT> output_bundle) {

    long offset = (in.ts - this->start).total_microseconds() \
        % (this->window_size).total_microseconds();

    output_bundle->add_value(
          Window(in.ts - microseconds(offset), this->window_size),
          text_to_score(in.data));

    return 1;
  }

private:

	long text_to_score(string_range const & str) {


		xzl_assert(str.data);
		xzl_assert(str.len < 1024); /* sane */
		string text (str.data, str.len);

#if 0 /* dbg only -- are we getting txt? */
		cout << text << endl;
#endif


		/* split @s into tokens and look up them in the keyword map.
		 * this is fast.
		 * cons: missing keywords.
		 * - different forms of a word? luckily, the input file already contain them:
		 * e.g. kill, killed, kills
		 * - multi words keyword. e.g. "cool stuff"
		 * XXX we may have to build a separate list for them and use each of them to
		 * search into the text.
		 */

		long score = 0;
		std::vector<std::string> tokens;
		boost::split(tokens, text, boost::is_any_of("\t "));
		xzl_assert(tokens.size() != 0);
		for (auto & token : tokens) {
			auto it = this->keywords.find(token);
			if (it != this->keywords.end()) {
				score += it->second;
			}
		}

#if 0		/* dbg -- can we produce scores? */
		if (score > 2) {
			cout << "score: "  << score << endl;
			cout << text << endl;
		}
#endif

		return score;
	}
};

#endif /* WINMAPPER_H_ */
