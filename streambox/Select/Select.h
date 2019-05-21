#ifndef SELECT_H
#define SELECT_H

// drop (internally mask) some records based some criteria
// @Input: element type.
template <typename InputT>
class Select : public PTransform {
public:
  Select(string name)
    : PTransform(name) { }

  PValue* apply(PValue* input) override {
    return PCollection::createPrimitiveOutputInternal(input->getPipeline());
  }

  // cannot do static virtual
  //  static virtual bool do_select(Record const & rec) = 0;
};

#endif /* SELECT_H */
