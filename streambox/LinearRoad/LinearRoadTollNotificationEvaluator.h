#ifndef LRB_TOLL_NOTIFICATION_EVAL_H
#define LRB_TOLL_NOTIFICATION_EVAL_H

#include <sstream>
#include "Values.h"
#include "core/SingleInputTransformEvaluator.h"
#include "LinearRoad/LinearRoadTollNotification.h"



template <typename InputT, typename OutputT,
  template<class> class BundleT>

class LinearRoadTollNotificationEvaluator
    : public SingleInputTransformEvaluator<
    LinearRoadTollNotification<InputT, OutputT, BundleT>, BundleT<InputT>,
    BundleT<OutputT>>
{

  using InputBundleT = BundleT<InputT>;
  using OutputBundleT = BundleT<OutputT>;
  using TransformT = LinearRoadTollNotification<InputT, OutputT, BundleT>;

public:
  bool evaluateSingleInput (TransformT* trans,
                            shared_ptr<InputBundleT> input_bundle,
                            shared_ptr<OutputBundleT> output_bundle) override {

	size_t pred1 = 0;
	size_t pred2 = 0;
	size_t pred3 = 0;
	size_t pred4 = 0;

	size_t toll_out_less_40 = 0;
	size_t toll_out_acc = 0;
	size_t toll_seq = 0;

    size_t toll_cnt = 0;
    size_t toll_calc = 0;
    unsigned long consumed = 0;
    for (auto && it = input_bundle->begin(); it != input_bundle->end(); ++it) {
      consumed ++;
      trans->record_counter_.fetch_add(1, std::memory_order_relaxed);

      auto rec = it->data;
      int current_pos = rec.m_iPos;
      int vid = rec.m_iVid;
      int time = rec.m_iTime;
      int speed = rec.m_iSpeed;
      int seq = rec.m_iPos/SEGMENT_LENGTH;

      int old_pos = -1;
      auto stop_map_it = trans->stop_map.find(vid);
      if (stop_map_it != trans->stop_map.end()) {
        old_pos = stop_map_it->second.pos;
      } else {
        trans->stop_map[vid] = StopMap();
      }

      bool possible_accident = false;
      /* for (auto map_it = trans->stop_map.begin(); map_it != trans->stop_map.end(); ++ map_it) { */
      /*   std::cout << "1: " << map_it->first << " -> " << map_it->second.toString() << std::endl; */
      /* } */
      if (speed == 0) {  /* car stops on the same position */
        if (old_pos == current_pos) {  /* stops at same position */
          trans->stop_map[vid].count ++;
          /* if 4 stops in row */
          if (trans->stop_map[vid].count == 4) {
            possible_accident = true;
          }
        }
      }
      // std::cout << "possible_accident: " << possible_accident << std::endl;

      if (old_pos != current_pos) {
        trans->stop_map[vid].pos = current_pos;
        trans->stop_map[vid].count = 1;
      }
      /* for (auto map_it = trans->stop_map.begin(); map_it != trans->stop_map.end(); ++ map_it) { */
      /*   std::cout << "2: " << map_it->first << " -> " << map_it->second.toString() << std::endl; */
      /* } */
      if (possible_accident) {  /* car has stopped 4 times in a row now. */
        trans->acc_dec->addStop(vid, current_pos, time);
      }
      // std::cout << "accident dectection size: " << trans->acc_dec->size() << std::endl;
      /* if car drives again, has speed. */
      if (speed != 0) {
        pred1 ++;
        if (trans->acc_dec->isInAccidentMapEntryForPos(old_pos)) {
          pred2 ++;
          if (trans->acc_dec->isVedAtPos(vid, old_pos)) {
            pred4 ++;
            trans->acc_dec->removeEntry(vid, old_pos);
          }
        }
      }

      trans->avg_seg_speed_map[seq].value += speed;
      trans->avg_seg_speed_map[seq].count ++;

      if (time % 300 == 0 && trans->time_of_last_toll != time && time != 0) {

        toll_cnt ++;
        for (size_t i = 0; i < SEGMENT_NUM; ++ i) {
          toll_seq ++;
          size_t current_avg_speed = trans->avg_seg_speed_map[i].value /
                                     trans->avg_seg_speed_map[i].count;

          // XXX: less than 40 ?
          if (current_avg_speed > 40) {
            toll_out_less_40 ++;
            bool accident = trans->acc_dec->isAccidentInSeq(seq);
            if (! accident) {
              toll_calc ++;
            } else {
              toll_seq ++;
            }
          } else {
          }
        }
        trans->time_of_last_toll = time;

        /*   /\* size_t avg_speed = trans->avg_speed_maps[it->data.m_iVid]->getAvgSpeed(); *\/ */
        /*   /\* OutputT out = OutputT(it->data.m_iType, *\/ */
        /*   /\*                       it->data.m_iTime, *\/ */
        /*   /\*                       it->data.m_iVid, *\/ */
        /*   /\*                       0, /\\* emit *\\/ *\/ */
        /*   /\*                       avg_speed, *\/ */
        /*   /\*                       0 *\/ */
        /*   /\*                       ); *\/ */

        /*   /\* auto out_record = Record<OutputT>(out, it->ts); *\/ */
        /*   /\* // auto out_record = Record<OutputT>((OutputT)it->data, it->ts); *\/ */
        /*   /\* std::cout << out_record.data.toString() << std::endl; *\/ */
        /*   /\* output_bundle->add_record(out_record); *\/ */
        /*   /\* // update timestamp of last toll *\/ */
        /*   /\* trans->time_of_last_toll = it->data.m_iTime; *\/ */
        /* } */
      }
    }

	/* stringstream ss; */
	/* ss << "processed " << consumed << " events and callAddStop=" << trans->acc_dec->addCallCnt << " calLDelCnt=" << trans->acc_dec->deleteCallCnt */
    /*    << " tollCnt=" << toll_cnt << " tollCalc=" << toll_calc << " addFirst=" << trans->acc_dec->addFirst << " addSecond=" << trans->acc_dec->addSecond  << " addNew=" << trans->acc_dec->addNew  << " addAlreadyFirst=" <<trans->acc_dec-> addAlreadyFirst */
    /*    << " addAlreadySecond=" << trans->acc_dec->addAlreadySecond  << " removeFirst=" << trans->acc_dec->removeFirst << " removeSecond=" << trans->acc_dec->removeSecond  << " erase=" << trans->acc_dec->erase */
    /*    << " couldNotBe=" << trans->acc_dec->couldNotBe */
    /*    << " pred1=" << pred1 << " pred2=" << pred2 << " pred4=" << pred4 */
    /*    << " tollOutLess40=" << toll_out_less_40 << " tollOutAcc=" << toll_out_acc << " tollSeq=" << toll_seq */
    /*    << endl; */
	/* cout << ss.str(); */

    // return false;
    return true;
  }

LinearRoadTollNotificationEvaluator(int node)
    : SingleInputTransformEvaluator<TransformT, InputBundleT, OutputBundleT>(node) { }
};


#endif // LRB_TOLL_NOTIFICATION_EVAL_H
