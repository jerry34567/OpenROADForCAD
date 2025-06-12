// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2025-2025, The OpenROAD Authors

#include "BaseMove.hh"
#include <climits>
#include <cmath>
#include <random>
namespace rsz {

using std::string;

using sta::ArcDelay;
using sta::Instance;
using sta::InstancePinIterator;
using sta::LoadPinIndexMap;
using sta::Net;
using sta::NetConnectedPinIterator;
using sta::Path;
using sta::PathExpanded;
using sta::Pin;
using sta::RiseFallBoth;
using sta::Slack;
using sta::Slew;
using sta::Vertex;

class SABufferMove : public BaseMove
{
 public:
  using BaseMove::BaseMove;
  const MinMax* min_ = MinMax::min();
  const MinMax* max_ = MinMax::max();

  bool doMove(const Path* drvr_path,
              int drvr_index,
              Slack drvr_slack,
              PathExpanded* expanded,
              float setup_slack_margin) override;
  bool myAddBuffer(const Path* drvr_path,
              int drvr_index,
              Slack drvr_slack,
              PathExpanded* expanded,
              float setup_slack_margin);

  const char* name() override { return "SABufferMove"; }

  // void rebufferNet(const Pin* drvr_pin);
 private:
  int rebuffer_net_count_ = 0;
  LibertyPort* drvr_port_ = nullptr;

  int rebuffer(const Pin* drvr_pin, int drvr_index);

  void annotateLoadSlacks(BufferedNetPtr& bnet, Vertex* root_vertex);
  BufferedNetPtr rebufferForTiming(const BufferedNetPtr& bnet, const Pin* drvr_pin);
  BufferedNetPtr recoverArea(const BufferedNetPtr& bnet,
                             sta::Delay slack_target,
                             float alpha);

  void debugCheckMultipleBuffers(Path* path, PathExpanded* expanded);
  bool hasTopLevelOutputPort(Net* net);

  int rebufferTopDown(const BufferedNetPtr& choice,
                      Net* net,
                      int level,
                      Instance* parent,
                      odb::dbITerm* mod_net_drvr,
                      odb::dbModNet* mod_net,
                      int drvr_index);
  BufferedNetPtr addWire(const BufferedNetPtr& p,
                         const Point& wire_end,
                         int wire_layer,
                         int level);
  void addBuffers(BufferedNetSeq& Z1,
                  int level,
                  bool area_oriented = false,
                  sta::Delay threshold = 0);
  float bufferInputCapacitance(LibertyCell* buffer_cell,
                               const DcalcAnalysisPt* dcalc_ap);
  Delay bufferDelay(LibertyCell* cell, const RiseFallBoth* rf, float load_cap);
  std::tuple<sta::Delay, sta::Delay> drvrPinTiming(const BufferedNetPtr& bnet);
  Slack slackAtDriverPin(const BufferedNetPtr& bnet);
  Slack slackAtDriverPin(const BufferedNetPtr& bnet, int index);

  Delay requiredDelay(const BufferedNetPtr& bnet);
  
  // Calculate total negative slack for a subtree
  Slack calculateSubtreeTns(BufferedNetPtr& bnet, const Pin* drvr_pin);

  // For rebuffering
  Path* arrival_paths_[RiseFall::index_count];
  void removeBuffer(Instance* buffer);
  double getT(){return _T;}
  double getr(){return _r;}
  double getcost(){return _cost;}
  double getcost_prev(){return _cost_prev;}
  double getcost_best(){return _cost_best;}
  size_t getiter_num(){return _iter_num;}
  void setT(double T){_T = T;}
  void setseed(size_t seed){srand(seed);}
  void setr(double r){_r = r;}
  void setIter_num(size_t iter_num){_iter_num = iter_num;}
  void setFrozen_point(double frozen_point){_frozen_point = frozen_point;}
  void setcost(double cost){_cost = cost;}
  void setcost_prev(double cost_prev){_cost_prev = cost_prev;}
  void setcost_best(double cost_best){_cost_best = cost_best;}
  bool update(double cost_new){
      _cost = cost_new;
      double delta_cost = -(cost_new - _cost_prev);
      if(delta_cost <= 0)
      {
          consecutive_reject_num = 0;
          // std::cout << "delta <= 0, accept, new_cost = " << cost_new << "old_cost = " << _cost_prev << std::endl;
          _cost_prev = cost_new;
          prev_buf.clear();
            for(auto it = current_buf.begin(); it != current_buf.end(); it++) prev_buf.push_back(!((*it).empty()));
          if(cost_new > _cost_best)
          {
            // std::cout << "best cost update = " << cost_new << std::endl;
            _cost_best = cost_new;
            best_buf.clear();
            for(auto it = current_buf.begin(); it != current_buf.end(); it++) best_buf.push_back(!((*it).empty()));
          }
          // for(auto jj : current_buf)
          //   std::cout << (!(jj.empty()));
          // std::cout << std::endl;
          return true;
      }
      else{
          double p = exp(-delta_cost*1000000000000/_T);
          // std::cout << "p = " << p << std::endl;
          bool ret = ((double)rand()/RAND_MAX < p);
          if(ret) // accept
          {
            consecutive_reject_num = 0;
            // std::cout << "delta > 0, accept, new_cost = " << cost_new << "old_cost = " << _cost_prev << std::endl;            _cost_prev = cost_new;
            prev_buf.clear();
            for(auto it = current_buf.begin(); it != current_buf.end(); it++) prev_buf.push_back(!((*it).empty()));
            if(cost_new > _cost_best)
            {
              _cost_best = cost_new;
              best_buf.clear();
              for(auto it = current_buf.begin(); it != current_buf.end(); it++) best_buf.push_back(!((*it).empty()));
            }
          }
          else // reject
          {
            ++consecutive_reject_num;
            // std::cout << "reject, bad_cost = " << cost_new << std::endl;
            // std::cout << "restoring..."<< std::endl;
            restore_buf(true);
            // std::cout << "restore_cost = " << _cost << std::endl;
          }
          // for(auto jj : current_buf)
          //   std::cout << (!(jj.empty()));
          // std::cout << std::endl;
          return ret;
      }
  }
  void cooling(){_T *= _r;}
  void perturb();
  bool isFrozen(){return _T < _frozen_point;}
  void restore_buf(bool is_prev);
  void reset();
  void myinit();
  double getTns(){
    resizer_->updateParasitics();
    sta_->findRequireds();
    return sta_->totalNegativeSlack(max_);
  }
  void final_result();
  double _T;
  double _frozen_point;
  double _r;
  double _cost_prev;
  double _cost_best;
  double _cost;
  PathExpanded* complete_path;
  size_t curr_buff_num;
  size_t _iter_num;
  vector<vector<Instance*>> current_buf; // {path_index, buffer_instance}
  vector<bool> best_buf; // == 1 -> has buffer
  vector<bool> prev_buf;
  vector<vector<Instance*>> orig_buf;
  size_t consecutive_reject_num;
  size_t max_consecutive_reject_num;
};

}  // namespace rsz
