#ifndef SIMPLEMAPPER_H
#define SIMPLEMAPPER_H

#include "Mapper/Mapper.h"
#include "Values.h"

//template <typename InputT, typename KVPair>
template <typename InputT = long, typename KVPair = pair<long, long>>
class SimpleMapper : public Mapper<InputT> {
public:
	SimpleMapper(string name = "dumb_mapper") : Mapper<InputT>(name) { }

	//static Record<KVPair> do_map(Record<InputT> const & in) {
	Record<KVPair> do_map(Record<InputT> const & in, int i) {
		//std::cout << "SimpleMapper do_map" << std::endl;
		return Record<KVPair>(
				//        KVPair(in.data / 100, in.data),  // an arbitrary mapper function
				//        KVPair(cnt, in.data),  // an arbitrary mapper function
				//KVPair(in.data, 2),
				KVPair(in.data, i),
				in.ts);
	}

	 
  	void ExecEvaluator(int nodeid, EvaluationBundleContext *c,
  		shared_ptr<BundleBase> bundle_ptr = nullptr) override;

};

#endif /* SIMPLEMAPPER_H */
