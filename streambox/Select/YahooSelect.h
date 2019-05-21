#ifndef YAHOOSELECT_H
#define YAHOOSELECT_H

// the "filtering" transform
template <typename InputT>
class YahooSelect : public Select<InputT> {
  using RecordT = Record<InputT>;
public:
  YahooSelect(string name = "yahoo_select") : Select<InputT>(name) { }

  // being static, there's not polymorphsim cost.
  inline static bool do_select(RecordT const & rec) {
    if (rec.ad_type.compare("view") != 0)
      return false;
    else
      return true;
  }
};

#endif /* YAHOOSELECT_H */
