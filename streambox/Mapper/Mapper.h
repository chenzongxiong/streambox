#ifndef MAPPER_H
#define MAPPER_H

#include "core/Transforms.h"

// input: a stream of records
// output: a stream of KV pairs (mostly for testing WinGBK)
// @InputT: the element type of input
template <typename InputT>
class Mapper: public PTransform {
	public:
		Mapper(string name) : PTransform(name) { }
};

/*
   template <typename InputT, typename KVPair>
   class SimpleMapper : public Mapper<InputT> {
   public:
   SimpleMapper(string name = "dumb_mapper") : Mapper<InputT>(name) { }

   static Record<KVPair> do_map(Record<InputT> const & in) {
   return Record<KVPair>(
//        KVPair(in.data / 100, in.data),  // an arbitrary mapper function
//        KVPair(cnt, in.data),  // an arbitrary mapper function
KVPair(in.data, 1),
in.ts);
}
};
 */
#endif // MAPPER_H
