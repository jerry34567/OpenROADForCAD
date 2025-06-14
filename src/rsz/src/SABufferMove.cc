// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2025-2025, The OpenROAD Authors

#include "SABufferMove.hh"

#include <algorithm>
#include <cmath>
#include <cstdio>

#include "BaseMove.hh"

namespace rsz {

using std::string;

using utl::RSZ;

using sta::ArcDelay;
using sta::Instance;
using sta::InstancePinIterator;
using sta::LoadPinIndexMap;
using sta::Net;
using sta::NetConnectedPinIterator;
using sta::Path;
using sta::PathExpanded;
using sta::Pin;
using sta::Slack;
using sta::Slew;
using sta::TimingArc;
using sta::Vertex;

bool SABufferMove::doMove(const Path* drvr_path,
                        int drvr_index,
                        Slack drvr_slack,
                        PathExpanded* expanded,
                        float setup_slack_margin)
{
    complete_path = expanded;
    curr_buff_num = 0;
    std::cout << "SA isUseFul!" << std::endl;
    myinit();
    while(!isFrozen())
    {
        for(int i = 0, ni = getiter_num(); i < ni; ++i)
        {
            perturb();
            update(getTns());
            if(consecutive_reject_num > max_consecutive_reject_num)
            {
                // std::cout << "too much rejection" << std::endl;
                setT(1);
            }
        }
        cooling();
        // std::cout << "cooling..." << std::endl;
    }
    // std::cout << "final result" << std::endl; 
    final_result();
    reset();
    return true;
}
bool SABufferMove::myAddBuffer(const Path* drvr_path,
                        int drvr_index,
                        Slack drvr_slack,
                        PathExpanded* expanded,
                        float setup_slack_margin)
{
  Vertex* drvr_vertex = drvr_path->vertex(sta_);
  const Pin* drvr_pin = drvr_vertex->pin();
  Instance* drvr_inst = network_->instance(drvr_pin);

  
  // Rebuffer blows up on large fanout nets.
//   if (fanout >= rebuffer_max_fanout_) {
//     return false;
//   }
  const bool tristate_drvr = resizer_->isTristateDriver(drvr_pin);
  if (tristate_drvr) {
    return false;
  }
  const Net* net = db_network_->dbToSta(db_network_->flatNet(drvr_pin));
  if (resizer_->dontTouch(net)) {
    return false;
  }
  dbNet* db_net = db_network_->staToDb(net);
  if (db_net->isConnectedByAbutment()) {
    return false;
  }

  const int rebuffer_count = rebuffer(drvr_pin, drvr_index);
  if (rebuffer_count > 0) {
    debugPrint(logger_,
               RSZ,
               "repair_setup",
               3,
               "rebuffer {} inserted {}",
               network_->pathName(drvr_pin),
               rebuffer_count);
    debugPrint(logger_,
               RSZ,
               "opt_moves",
               1,
               "ACCEPT buffer {} inserted {}",
               network_->pathName(drvr_pin),
               rebuffer_count);
    addMove(drvr_inst, rebuffer_count);
  } else {
    debugPrint(logger_,
               RSZ,
               "opt_moves",
               3,
               "REJECT buffer {} inserted {}",
               network_->pathName(drvr_pin),
               rebuffer_count);
  }
  return rebuffer_count > 0;
}

void SABufferMove::debugCheckMultipleBuffers(Path* path, PathExpanded* expanded)
{
  if (expanded->size() > 1) {
    const int path_length = expanded->size();
    const int start_index = expanded->startIndex();
    for (int i = start_index; i < path_length; i++) {
      const Path* path = expanded->path(i);
      const Pin* path_pin = path->pin(sta_);
      if (i > 0 && network_->isDriver(path_pin)
          && !network_->isTopLevelPort(path_pin)) {
        const TimingArc* prev_arc = path->prevArc(sta_);
        printf("repair_setup %s: %s ---> %s \n",
               prev_arc->from()->libertyCell()->name(),
               prev_arc->from()->name(),
               prev_arc->to()->name());
      }
    }
  }
  printf("done\n");
}
void SABufferMove::removeBuffer(Instance* buffer)
{
  debugPrint(logger_,
             RSZ,
             "repair_setup",
             3,
             "remove_buffer{}",
             network_->pathName(buffer));
  addMove(buffer);

  LibertyCell* lib_cell = network_->libertyCell(buffer);
  LibertyPort *in_port, *out_port;
  lib_cell->bufferPorts(in_port, out_port);

  Pin* in_pin = db_network_->findPin(buffer, in_port);
  Pin* out_pin = db_network_->findPin(buffer, out_port);

  // Hierarchical net handling
  odb::dbModNet* op_modnet = db_network_->hierNet(out_pin);

  odb::dbNet* in_db_net = db_network_->flatNet(in_pin);
  odb::dbNet* out_db_net = db_network_->flatNet(out_pin);
  if (in_db_net == nullptr || out_db_net == nullptr) {
    return;
  }
  // in_net and out_net are flat nets.
  Net* in_net = db_network_->dbToSta(in_db_net);
  Net* out_net = db_network_->dbToSta(out_db_net);

  bool out_net_ports = hasPort(out_net);
  Net *survivor, *removed;
  if (out_net_ports) {
    survivor = out_net;
    removed = in_net;
  } else {
    // default or out_net_ports
    // Default to in_net surviving so drivers (cached in dbNetwork)
    // do not change.
    survivor = in_net;
    removed = out_net;
  }
  debugPrint(
      logger_, RSZ, "remove_buffer", 1, "remove {}", db_network_->name(buffer));

  odb::dbNet* db_survivor = db_network_->staToDb(survivor);
  odb::dbNet* db_removed = db_network_->staToDb(removed);
  if (db_removed) {
    db_survivor->mergeNet(db_removed);
  }
  sta_->disconnectPin(in_pin);
  sta_->disconnectPin(out_pin);
  sta_->deleteInstance(buffer);
  if (removed) {
    sta_->deleteNet(removed);
  }

  // Hierarchical case supported:
  // moving an output hierarchical net to the input pin driver.
  // During canBufferRemove check (see above) we require that the
  // input pin driver is in the same module scope as the output hierarchical
  // driver
  //
  if (op_modnet) {
    debugPrint(logger_,
               RSZ,
               "remove_buffer",
               1,
               "Handling hierarchical net {}",
               op_modnet->getName());
    Pin* driver_pin = nullptr;
    db_network_->getNetDriverParentModule(in_net, driver_pin, true);
    db_network_->connectPin(driver_pin, db_network_->dbToSta(op_modnet));
  }

  resizer_->parasitics_invalid_.erase(removed);
  resizer_->parasiticsInvalid(survivor);
  resizer_->updateParasitics();
}

void SABufferMove::restore_buf(bool is_prev)
{
    // curr_buff_num = 0;
    vector<bool>& new_buf = is_prev ? prev_buf : best_buf;
    for(int i = 0, ni = complete_path->size(); i < ni; i++)
    {
        if(current_buf.at(i).empty() == new_buf.at(i))
        {
            if(!new_buf.at(i))
            {
                for(auto buf_inst: current_buf.at(i))
                    removeBuffer(buf_inst);
                current_buf.at(i).clear();
            }
            else
            {
                myAddBuffer(complete_path->path(i), i, sta_->vertexSlack(complete_path->path(i)->vertex(sta_), max_), complete_path, 0);
            }
        }
    }
    
}
void SABufferMove::perturb()
{
    int fanout = 0;
    int toggle_path = 0;
    while(fanout < 1)
    {
        toggle_path = rand() % complete_path->size();
        Vertex* drvr_vertex = complete_path->path(toggle_path)->vertex(sta_);
        fanout = this->fanout(drvr_vertex);
    }
    // std::cout << "perturbed path: " << toggle_path << std::endl;
    if(current_buf.at(toggle_path).empty())
    {
        // std::cout << "add Buffer" << std::endl;
        myAddBuffer(complete_path->path(toggle_path), toggle_path, sta_->vertexSlack(complete_path->path(toggle_path)->vertex(sta_), max_), complete_path, 0);
        // curr_buff_num++;
    }
    else
    {
        // std::cout << "remove Buffer" << std::endl;
        for(auto buf_inst: current_buf.at(toggle_path))
            removeBuffer(buf_inst);
        current_buf.at(toggle_path).clear();
        // curr_buff_num--;
    }
    return;
}
void SABufferMove::myinit()
{
    setIter_num(20);
    setT(1000);
    setr(0.8);
    setFrozen_point(1);
    consecutive_reject_num = 0;
    max_consecutive_reject_num = 10;
    double curr_tns = getTns();
    _cost = curr_tns;
    _cost_prev = curr_tns;
    _cost_best = curr_tns;
    current_buf.clear();
    best_buf.clear();
    prev_buf.clear();
    current_buf.resize(complete_path->size());
    best_buf.resize(complete_path->size());
    prev_buf.resize(complete_path->size());
    std::fill(current_buf.begin(), current_buf.end(), vector<Instance*>());
    std::fill(best_buf.begin(), best_buf.end(), false);
    std::fill(prev_buf.begin(), prev_buf.end(), false);
    return;
}
void SABufferMove::final_result()
{
    // remember orignal
    // for(int i = 0, ni = complete_path->size(); i < ni; i++)
    // {
    //     if(current_buf.at(i).empty() == best_buf.at(i))
    //     {
    //         if(!best_buf.at(i))
    //         {
    //             for(auto buf_inst: current_buf.at(i))
    //                 removeBuffer(buf_inst);
    //         }
    //         else
    //         {
    //             myAddBuffer(complete_path->path(i), i, sta_->vertexSlack(complete_path->path(i)->vertex(sta_), max_), complete_path, 0);
    //         }
    //     }
    // }
    restore_buf(false);
    return;
}
void SABufferMove::reset()
{
    current_buf.clear();
    best_buf.clear();
    prev_buf.clear();
    orig_buf.clear();
}
}  // namespace rsz
