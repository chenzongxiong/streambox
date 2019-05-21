// @T: element type, e.g. long
template<class T>
class UnboundedPseudoEvaluator : public TransformEvaulator<UnboundedPseudo<T>> {
public:
  /* an infi loop. emit a bundle periodically.
   * note: since IO bound, this cannot be executed as TBB or Cilk's task
   * which will block the underlying worker thread. */
  void evaluate(UnboundedPseudo<T>* t, EvaluationBundleContext* c) override {

    boost::posix_time::time_duration delta = milliseconds(t->interval_ms);

    auto out = t->getFirstOutput();
    auto out2 = t->getSecondOutput();
    assert(out);

    while (true) {

      W(" ----------- to sleep %d ms", t->interval_ms);

      usleep(t->interval_ms * 1000);

      // the # of bundles to emit before pause and advance watermark
//      const int bundle_count = 2 * std::thread::hardware_concurrency();
      const int bundle_count = 2 * c->num_workers;

      vector<std::future<int>> fs;

      // issue multi bundles at a time. each bundle spawns a new task
      for (int i = 0; i < bundle_count; i++) {
        auto f = c->tp.push([this, i, delta, out, c](int id) { // lambda
          shared_ptr<RecordBitmapBundle<T>>
            bundle(make_shared<RecordBitmapBundle<T>>(max_bundle_size));
          assert(bundle);

          /* all records in the bundle have same ts */
          for (int j = 0; j < max_bundle_size; j ++) {
              bundle->add_record(Record<T>(genRandomElement(),
                  current_ts + delta * i));
          }

          out->depositOneBundle(bundle);
          c->SpawnConsumer(out);

          return 1;
        }); // end of lambda
        fs.push_back(std::move(f));  // future has to be moved.

        /* a secondary output stream */
        if (out2) {
            auto f = c->tp.push([this, i, delta, out2, c](int id) { // lambda
              shared_ptr<RecordBitmapBundle<T>>
                bundle(make_shared<RecordBitmapBundle<T>>(max_bundle_size));
              assert(bundle);

              for (int j = 0; j < max_bundle_size; j ++) {
                  bundle->add_record(Record<T>(genRandomElement1(),
                      current_ts + delta * i));
              }

              out2->depositOneBundle(bundle);
              c->SpawnConsumer(out2);

              return 1;
            }); // end of lambda
            fs.push_back(std::move(f));  // future has to be moved.
        }
      }

      // make sure all data have left the current evaluator before
      // advancing the watermark
      for (auto&& fu : fs) {
          VV("fu get...");
          fu.get();
      }

      // expensive?
      current_ts = current_ts + delta * bundle_count;

//      current_ts = current_ts + milliseconds(5000);
      current_ts = current_ts + milliseconds(t->session_gap_ms);

//      c->SeedWaterMark(current_ts - seconds(60)); // lag
      c->SeedWaterMark(current_ts);
    }
  }

  ~UnboundedPseudoEvaluator() { }

//  UnboundedPseudoEvaluator(int max_bundle_size = 1024000)
  UnboundedPseudoEvaluator(int max_bundle_size = 1024)
  : max_bundle_size (max_bundle_size) { }

private:
  int max_bundle_size;
  ptime current_ts = Window::epoch; // easy to debug
//  atomic<long> _cnt;

  inline T genRandomElement() { assert(0); } // generic, not impl
  inline T genRandomElement1() { assert(0); } // generic, not impl
};

/* -- template specialization -- */
template<>
inline long UnboundedPseudoEvaluator<long>::genRandomElement() {
//  return rand(); // xzl: can be expensive.
//  return 1234;

#if 0
  static atomic<long> cnt (0);
  long reset = 0, m = 40;
  if (cnt.fetch_add(1) == m) {
      cnt.compare_exchange_strong(m, reset);
  }
  return cnt;
#endif

  static atomic<long> cnt (0);
  return cnt.fetch_add(1);
}

template<>
inline long UnboundedPseudoEvaluator<long>::genRandomElement1() {
  static atomic<long> cnt (0);
  return cnt.fetch_add(1);
}


template<>
inline tuple<long, float> UnboundedPseudoEvaluator<tuple<long, float>>::genRandomElement() {
  return make_tuple(rand(), (float) rand() / RAND_MAX);
}
