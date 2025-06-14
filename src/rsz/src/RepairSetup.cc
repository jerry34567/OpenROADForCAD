// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2019-2025, The OpenROAD Authors

#include "RepairSetup.hh"

#include <algorithm>
#include <cfloat>
#include <cmath>
#include <cstddef>
#include <cstring>
#include <limits>
#include <memory>
#include <optional>
#include <random>
#include <sstream>
#include <string>

#include "BaseMove.hh"
#include "BufferMove.hh"
#include "CloneMove.hh"
#include "GAMove.hh"
#include "SizeDownMove.hh"
#include "SizeUpMove.hh"
#include "SplitLoadMove.hh"
#include "SwapPinsMove.hh"
#include "UnbufferMove.hh"
#include "rsz/Resizer.hh"
#include "sta/Corner.hh"
#include "sta/DcalcAnalysisPt.hh"
#include "sta/Fuzzy.hh"
#include "sta/Graph.hh"
#include "sta/GraphDelayCalc.hh"
#include "sta/InputDrive.hh"
#include "sta/Liberty.hh"
#include "sta/Parasitics.hh"
#include "sta/PathExpanded.hh"
#include "sta/PortDirection.hh"
#include "sta/Sdc.hh"
#include "sta/TimingArc.hh"
#include "sta/Units.hh"
#include "sta/VerilogWriter.hh"
#include "utl/Logger.h"
#include "sta/Path.hh"

namespace rsz {

using std::max;
using std::pair;
using std::string;
using std::vector;
using utl::RSZ;

using sta::Edge;
using sta::fuzzyEqual;
using sta::fuzzyGreater;
using sta::fuzzyGreaterEqual;
using sta::fuzzyLess;
using sta::GraphDelayCalc;
using sta::InstancePinIterator;
using sta::NetConnectedPinIterator;
using sta::PathExpanded;
using sta::PathSeq;
using sta::Slew;
using sta::VertexOutEdgeIterator;
using sta::VertexPathIterator;

RepairSetup::RepairSetup(Resizer* resizer) : resizer_(resizer)
{
}

void RepairSetup::init()
{
  logger_ = resizer_->logger_;
  dbStaState::init(resizer_->sta_);
  db_network_ = resizer_->db_network_;

  initial_design_area_ = resizer_->computeDesignArea();
  
  // Clear TNS window for fresh start
  while (!shuffle_tns_windows_.empty()) {
    shuffle_tns_windows_.pop();
  }
}

bool RepairSetup::repairSetup(const float setup_slack_margin,
                              const double repair_tns_end_percent,
                              const int max_passes,
                              const int max_repairs_per_pass,
                              const bool verbose,
                              const std::vector<MoveType>& sequence,
                              const bool skip_pin_swap,
                              const bool skip_gate_cloning,
                              const bool skip_gate_sizing,
                              const bool skip_size_down,
                              const bool skip_buffering,
                              const bool skip_split_load,
                              const bool skip_buffer_removal,
                              const bool skip_last_gasp,
                              const bool ga_enabled,
                              const bool shuffle_enabled)
{
  bool repaired = false;
  init();
  if (shuffle_enabled) {  
    initializeGaGates();
  }
  constexpr int digits = 3;
  max_repairs_per_pass_ = max_repairs_per_pass;
  resizer_->buffer_moved_into_core_ = false;
  ga_enabled_ = ga_enabled;
  shuffle_enabled_ = shuffle_enabled;
  if (!sequence.empty()) {
    move_sequence.clear();
    for (MoveType move : sequence) {
      switch (move) {
        case MoveType::BUFFER:
          if (!skip_buffering) {
            move_sequence.push_back(resizer_->buffer_move);
          }
          break;
        case MoveType::UNBUFFER:
          if (!skip_buffer_removal) {
            move_sequence.push_back(resizer_->unbuffer_move);
          }
          break;
        case MoveType::SWAP:
          if (!skip_pin_swap) {
            move_sequence.push_back(resizer_->swap_pins_move);
          }
          break;
        case MoveType::SIZE:
          if (!skip_size_down) {
            move_sequence.push_back(resizer_->size_down_move);
          }
          move_sequence.push_back(resizer_->size_up_move);
          break;
        case MoveType::SIZEUP:
          move_sequence.push_back(resizer_->size_up_move);
          break;
        case MoveType::SIZEDOWN:
          if (!skip_size_down) {
            move_sequence.push_back(resizer_->size_down_move);
          }
          break;
        case MoveType::CLONE:
          if (!skip_gate_cloning) {
            move_sequence.push_back(resizer_->clone_move);
          }
          break;
        case MoveType::SPLIT:
          if (!skip_split_load) {
            move_sequence.push_back(resizer_->split_load_move);
          }
          break;
      }
    }

  } else {
    move_sequence.clear();
    if (!skip_buffer_removal) {
      move_sequence.push_back(resizer_->unbuffer_move);
    }
    // TODO: Add size_down_move to the sequence if we want to allow
    // Always  have sizing
    if (!skip_gate_sizing) {
        move_sequence.push_back(resizer_->size_up_move);
    }
    if (!skip_pin_swap) {
      move_sequence.push_back(resizer_->swap_pins_move);
    }
    if (!skip_buffering) {
      move_sequence.push_back(resizer_->buffer_move);
    }
    if (!skip_gate_cloning) {
      move_sequence.push_back(resizer_->clone_move);
    }
    if (!skip_split_load) {
      move_sequence.push_back(resizer_->split_load_move);
    }
    if (ga_enabled) {
      move_sequence.push_back(resizer_->ga_move);
    }
  }

  string repair_moves = "Repair move sequence: ";
  for (auto move : move_sequence) {
    repair_moves += move->name() + string(" ");
  }
  logger_->info(RSZ, 100, repair_moves);

  // Sort failing endpoints by slack.
  const VertexSet* endpoints = sta_->endpoints();
//   vector<pair<Vertex*, Slack>> violating_ends;
  // logger_->setDebugLevel(RSZ, "repair_setup", 2);
  // Should check here whether we can figure out the clock domain for each
  // vertex. This may be the place where we can do some round robin fun to
  // individually control each clock domain instead of just fixating on fixing
  // one.
//   for (Vertex* end : *endpoints) {
//     const Slack end_slack = sta_->vertexSlack(end, max_);
//     if (end_slack < setup_slack_margin) {
//       violating_ends.emplace_back(end, end_slack);
//     }
//   }
//   std::stable_sort(violating_ends.begin(),
//                    violating_ends.end(),
//                    [](const auto& end_slack1, const auto& end_slack2) {
//                      return end_slack1.second < end_slack2.second;
//                    });
  PathSeq violating_paths;
  for (Vertex* end : *endpoints) {
    const Slack end_slack = sta_->vertexSlack(end, max_);
    VertexPathIterator path_iter(end, this);
    PathSeq paths;
    while (path_iter.hasNext()) {
      Path* path = path_iter.next();
      paths.push_back(path);
    }
    std::stable_sort(paths.begin(),
                     paths.end(),
                     [this](const auto& path1, const auto& path2) {
                       return path1->slack(sta_) < path2->slack(sta_);
                     });
    for (int i = 0; i < 5 && i < paths.size(); i++) {
      violating_paths.push_back(paths[i]);
    }
  }
  std::stable_sort(violating_paths.begin(),
                   violating_paths.end(),
                   [this](const auto& path1, const auto& path2) {
                     return path1->slack(sta_) < path2->slack(sta_);
                   });
  debugPrint(logger_,
             RSZ,
             "repair_setup",
             1,
             "Violating endpoints {}/{} {}%",
             violating_paths.size(),
             endpoints->size(),
             int(violating_paths.size() / double(endpoints->size()) * 100));

  if (!violating_paths.empty()) {
    logger_->info(RSZ,
                  94,
                  "Found {} endpoints with setup violations.",
                  violating_paths.size());
  } else {
    // nothing to repair
    logger_->metric("design__instance__count__setup_buffer", 0);
    logger_->info(RSZ, 98, "No setup violations found");
    return false;
  }

  int end_index = 0;
  int max_end_count = violating_paths.size() * repair_tns_end_percent;
  float initial_tns = sta_->totalNegativeSlack(max_);
  best_tns_ = initial_tns;
  float prev_tns = initial_tns;
  int num_viols = violating_paths.size();
  // Always repair the worst endpoint, even if tns percent is zero.
  max_end_count = max(max_end_count, 1);
  logger_->info(RSZ,
                99,
                "Repairing {} out of {} ({:0.2f}%) violating endpoints...",
                max_end_count,
                violating_paths.size(),
                repair_tns_end_percent * 100.0);

  // Ensure that max cap and max fanout violations don't get worse
  sta_->checkCapacitanceLimitPreamble();
  sta_->checkFanoutLimitPreamble();

  IncrementalParasiticsGuard guard(resizer_);
  int opto_iteration = 0;
  bool prev_termination = false;
  bool two_cons_terminations = false;
  printProgress(opto_iteration, false, false, false, num_viols);
  float fix_rate_threshold = inc_fix_rate_threshold_;
  if (!violating_paths.empty()) {
    min_viol_ = -violating_paths.back()->slack(sta_);
    max_viol_ = -violating_paths.front()->slack(sta_);
  }
  for (int j = 0; j < shuffle_times_; j++) {
    end_index = 0;
    while (!violating_paths.empty()) {
      auto it = violating_paths.begin();
      while (it != violating_paths.end()) {
        if ((*it)->slack(sta_) >= 0) {
          it = violating_paths.erase(it);
        } else {
          ++it;
        }
      }
      std::stable_sort(violating_paths.begin(),
                    violating_paths.end(),
                    [this](const auto& path1, const auto& path2) {
                      return path1->slack(sta_) < path2->slack(sta_);
                    });
      Path* end_path = violating_paths.front();
      violating_paths.erase(violating_paths.begin());
      std::cout << "path: " << end_path->vertex(sta_)->name(network_) << " " << end_path->slack(sta_) << std::endl;

      fallback_ = false;
      Vertex* end = end_path->vertex(sta_);
      Slack end_slack = sta_->vertexSlack(end, max_);
      Slack worst_slack;
      Vertex* worst_vertex;
      sta_->worstSlack(max_, worst_slack, worst_vertex);
      debugPrint(logger_,
                RSZ,
                "repair_setup",
                1,
                "{} slack = {} worst_slack = {}",
                end->name(network_),
                delayAsString(end_slack, sta_, digits),
                delayAsString(worst_slack, sta_, digits));
      end_index++;
      debugPrint(logger_,
                RSZ,
                "repair_setup",
                1,
                "Doing {} /{}",
                end_index,
                max_end_count);
      if (end_index > max_end_count) {
        // clang-format off
        debugPrint(logger_, RSZ, "repair_setup", 1, "{} end_index {} is larger than"
                  " max_end_count {}", end->name(network_), end_index,
                  max_end_count);
        // clang-format on
        break;
      }
      Slack prev_end_slack = end_slack;
      Slack prev_worst_slack = worst_slack;
      Delay prev_tns = sta_->totalNegativeSlack(max_);
      
      int pass = 1;
      int decreasing_slack_passes = 0;
      resizer_->journalBegin();
      while (pass <= max_passes) {
        opto_iteration++;
        if (verbose || opto_iteration == 1) {
          printProgress(opto_iteration, true, false, false, num_viols);
        }
        // if (terminateProgress(opto_iteration,
        //                       initial_tns,
        //                       prev_tns,
        //                       fix_rate_threshold,
        //                       end_index,
        //                       max_end_count)) {
        //   if (prev_termination) {
        //     // Abort entire fixing if no progress for 200 iterations
        //     two_cons_terminations = true;
        //   } else {
        //     prev_termination = true;
        //   }

        // Restore to previous good checkpoint
        // debugPrint(logger_,
        //            RSZ,
        //            "repair_setup",
        //            2,
        //            "Restoring best slack end slack {} worst slack {}",
        //            delayAsString(prev_end_slack, sta_, digits),
        //            delayAsString(prev_worst_slack, sta_, digits));
        // resizer_->journalRestore();
        // break;
      // }
        if (opto_iteration % opto_small_interval_ == 0) {
          prev_termination = false;
        }

        if (end_slack > setup_slack_margin) {
          --num_viols;
          if (pass != 1) {
            debugPrint(logger_,
                      RSZ,
                      "repair_setup",
                      2,
                      "Restoring best slack end slack {} worst slack {}",
                      delayAsString(prev_end_slack, sta_, digits),
                      delayAsString(prev_worst_slack, sta_, digits));
            resizer_->journalRestore();
          } else {
            resizer_->journalEnd();
          }
          // clang-format off
          debugPrint(logger_, RSZ, "repair_setup", 1, "bailing out at {}/{} "
                    "end_slack {} is larger than setup_slack_margin {}",
                    end_index, max_end_count, end_slack, setup_slack_margin);
          // clang-format on
          break;
        }
        // Path* end_path = sta_->vertexWorstSlackPath(end, max_);

        // bool changed = false;
        // VertexPathIterator* path_iter = sta_->vertexPathIterator(end, nullptr, max_);
        // while (path_iter->hasNext()) {
        //   Path *path = path_iter->next();
        //   Slack slack = path->slack(sta_);
        //   if (slack < setup_slack_margin) {
        //     bool b = repairPath(path, slack, setup_slack_margin);
        //     changed = changed || b;
        //   }
        // }
        const bool changed = repairPath(end_path, end_slack, setup_slack_margin);
        if (!changed) {
          if (pass != 1) {
            debugPrint(logger_,
                      RSZ,
                      "repair_setup",
                      2,
                      "No change after {} decreasing slack passes.",
                      decreasing_slack_passes);
            debugPrint(logger_,
                      RSZ,
                      "repair_setup",
                      2,
                      "Restoring best slack end slack {} worst slack {}",
                      delayAsString(prev_end_slack, sta_, digits),
                      delayAsString(prev_worst_slack, sta_, digits));
            resizer_->journalRestore();
          } else {
            resizer_->journalEnd();
          }
          // clang-format off
          debugPrint(logger_, RSZ, "repair_setup", 1, "bailing out {} no changes"
                    " after {} decreasing passes", end->name(network_),
                    decreasing_slack_passes);
          // clang-format on
          break;
        }
        resizer_->updateParasitics();
        sta_->findRequireds();
        end_slack = sta_->vertexSlack(end, max_);
        Delay tns = sta_->totalNegativeSlack(max_);
        sta_->worstSlack(max_, worst_slack, worst_vertex);
        const bool better = (fuzzyGreater(tns, prev_tns)
                            || (end_index != 1 && fuzzyEqual(tns, prev_tns)
                                && fuzzyGreater(end_slack, prev_end_slack)));
        debugPrint(logger_,
                  RSZ,
                  "repair_setup",
                  2,
                  "pass {} slack = {} worst_slack = {} tns {} {}",
                  pass,
                  delayAsString(end_slack, sta_, digits),
                  delayAsString(worst_slack, sta_, digits),
                  delayAsString(tns, sta_, digits),
                  better ? "save" : "");
        if (better) {
          if (end_slack > setup_slack_margin) {
            --num_viols;
          }
          prev_tns = tns;
          prev_end_slack = end_slack;
          prev_worst_slack = worst_slack;
          decreasing_slack_passes = 0;
          resizer_->journalEnd();
          // Progress, Save checkpoint so we can back up to here.
          resizer_->journalBegin();
          
        } else {
          fallback_ = true;
          // Allow slack to increase to get out of local minima.
          // Do not update prev_end_slack so it saves the high water mark.
          decreasing_slack_passes++;
          if (decreasing_slack_passes > decreasing_slack_max_passes_) {
            // Undo changes that reduced slack.
            debugPrint(logger_,
                      RSZ,
                      "repair_setup",
                      2,
                      "decreasing slack for {} passes.",
                      decreasing_slack_passes);
            debugPrint(logger_,
                      RSZ,
                      "repair_setup",
                      2,
                      "Restoring best end slack {} worst slack {}",
                      delayAsString(prev_end_slack, sta_, digits),
                      delayAsString(prev_worst_slack, sta_, digits));
            resizer_->journalRestore();
            // clang-format off
            debugPrint(logger_, RSZ, "repair_setup", 1, "bailing out {} decreasing"
                      " passes {} > decreasig pass limit {}", end->name(network_),
                      decreasing_slack_passes, decreasing_slack_max_passes_);
            // clang-format on
            break;
          }
        }

        // if (resizer_->overMaxArea()) {
        //   // clang-format off
        //   debugPrint(logger_, RSZ, "repair_setup", 1, "bailing out {} resizer"
        //              " over max area", end->name(network_));
        //   // clang-format on
        //   resizer_->journalEnd();
        //   break;
        // }
        if (end_index == 1) {
          end = worst_vertex;
        }
        pass++;
      }  // while pass <= max_passes
      if (verbose || opto_iteration == 1) {
        printProgress(opto_iteration, true, false, false, num_viols);
      }

      if (shuffle_enabled_) {
        const Slack current_tns = sta_->totalNegativeSlack(max_);
        if (best_tns_ < current_tns) {
          best_tns_ = current_tns;
          debugPrint(logger_,
                      RSZ,
                      "shuffle",
                      2,
                      "store best cells for shuffle");
          for (int i = 0; i < original_cells_.size(); i++) {
            original_cells_[i] = network_->libertyCell(ga_gates_[i].instance);
          }
        }
      }

      // Update TNS window for shuffle decision
      if (shuffle_enabled_) {
        // shuffle_tns_windows_.push(current_tns);
        
        // // Keep only last 10 TNS values
        // if (shuffle_tns_windows_.size() > 20) {
        //   shuffle_tns_windows_.pop();
        // }
        
        // // Check if we should shuffle based on TNS slope
        // if (shuffle_tns_windows_.size() == 20) {
          // Calculate slope between first and last TNS value
          // const Delay first_tns = shuffle_tns_windows_.front();
          // const Delay last_tns = shuffle_tns_windows_.back();
          
          // // Debug output to check values
          // debugPrint(logger_,
          //           RSZ,
          //           "shuffle",
          //           1,
          //           "Debug: window size={}, first_tns={}, last_tns={}, diff={}",
          //           shuffle_tns_windows_.size(),
          //           delayAsString(first_tns, sta_, 6),
          //           delayAsString(last_tns, sta_, 6),
          //           last_tns - first_tns);
          
          // const float tns_slope = static_cast<float>(first_tns - last_tns) / initial_tns;
          // std::cout << "tns_slope: " << tns_slope << std::endl;
          
          // Trigger shuffle if slope is too flat (not improving enough)
          // if (tns_slope < shuffle_tns_slope_) {
          //   restoreOriginalSizes();
          //   // Clear window for next cycle
          //   while (!shuffle_tns_windows_.empty()) {
          //     shuffle_tns_windows_.pop();
          //   }
            
            // Perform shuffle
            // for (int i = 0; i < original_cells_.size() * shuffle_percent_; i++) {
            //   int random_index = ga_random_generator_() % original_cells_.size();
            //   int max_size_index = ga_gates_[random_index].available_cells.size() - 1;
            //   std::uniform_int_distribution<int> size_dist(ga_gates_[random_index].current_size_index, max_size_index);
            //   int new_size_index = size_dist(ga_random_generator_);
            //   Instance* inst = ga_gates_[random_index].instance;
            //   LibertyCell* new_cell = ga_gates_[random_index].available_cells[new_size_index];
            //   resizer_->replaceCell(inst, new_cell, true);
            // }
            // resizer_->updateParasitics();
            // sta_->findRequireds();
            // Slack wns;
            // Vertex* worst_vertex;
            // // useful for CAD
            // sta_->worstSlack(max_, wns, worst_vertex);
            // const Slack tns_after_shuffle = sta_->totalNegativeSlack(max_);
            // debugPrint(logger_,
            //           RSZ,
            //           "shuffle",
            //           2,
            //           "Shuffled {} gates. Slope: {:.6f}, first: {}, last: {}, WNS: {} TNS: {}",
            //           static_cast<int>(original_cells_.size() * shuffle_percent_),
            //           tns_slope,
            //           delayAsString(first_tns, sta_, digits),
            //           delayAsString(last_tns, sta_, digits),
            //           delayAsString(wns, sta_, digits),
            //           delayAsString(tns_after_shuffle, sta_, digits));
            // const VertexSet* endpoints = sta_->endpoints();
            // violating_paths.clear();
            // for (Vertex* end : *endpoints) {
            //   const Slack end_slack = sta_->vertexSlack(end, max_);
            //   VertexPathIterator path_iter(end, this);
            //   // std::cout << "end: " << end->name(network_) << " " << end_slack << std::endl;
            //   while (path_iter.hasNext()) {
            //     Path* path = path_iter.next();
            //     // std::cout << "path: " << path->vertex(sta_)->name(network_) << " " << path->slack(sta_) << std::endl;
            //     if (path->slack(sta_) < setup_slack_margin) {
            //       violating_paths.push_back(path);
            //     }
            //   }
            // }
            // shuffle_tns_slope_ *= 0.8;
            // shuffle_percent_ *= 0.8;
        //   }
        // }
      }
      // if (two_cons_terminations) {
      //   // clang-format off
      //   debugPrint(logger_, RSZ, "repair_setup", 1, "bailing out of setup fixing"
      //              "due to no TNS progress for two opto cycles");
      //   // clang-format on
      //   break;
      // }
    }  // for each violating endpoint

    if (shuffle_enabled_ && j < shuffle_times_ - 1) {
      for (int i = 0; i < original_cells_.size() * shuffle_percent_; i++) {
        int random_index = ga_random_generator_() % original_cells_.size();
        int max_size_index = ga_gates_[random_index].available_cells.size() - 1;
        std::uniform_int_distribution<int> size_dist(ga_gates_[random_index].current_size_index, max_size_index);
        int new_size_index = size_dist(ga_random_generator_);
        Instance* inst = ga_gates_[random_index].instance;
        LibertyCell* new_cell = ga_gates_[random_index].available_cells[new_size_index];
        resizer_->replaceCell(inst, new_cell, true);
      }
      resizer_->updateParasitics();
      sta_->findRequireds();
      Slack wns;
      Vertex* worst_vertex;
      // useful for CAD
      sta_->worstSlack(max_, wns, worst_vertex);
      const Slack tns_after_shuffle = sta_->totalNegativeSlack(max_);
      debugPrint(logger_,
                RSZ,
                "shuffle",
                2,
                "Shuffled {} gates. WNS: {} TNS: {}",
                static_cast<int>(original_cells_.size() * shuffle_percent_),
                delayAsString(wns, sta_, digits),
                delayAsString(tns_after_shuffle, sta_, digits));
      const VertexSet* endpoints = sta_->endpoints();
      violating_paths.clear();
      for (Vertex* end : *endpoints) {
        const Slack end_slack = sta_->vertexSlack(end, max_);
        VertexPathIterator path_iter(end, this);
        // std::cout << "end: " << end->name(network_) << " " << end_slack << std::endl;
        while (path_iter.hasNext()) {
          Path* path = path_iter.next();
          // std::cout << "path: " << path->vertex(sta_)->name(network_) << " " << path->slack(sta_) << std::endl;
          if (path->slack(sta_) < setup_slack_margin) {
            violating_paths.push_back(path);
          }
        }
      }
    }
  }

  if (!skip_last_gasp) {
    // do some last gasp setup fixing before we give up
    OptoParams params(setup_slack_margin, verbose);
    params.iteration = opto_iteration;
    params.initial_tns = initial_tns;
    repairSetupLastGasp(params, num_viols);
  }

  printProgress(opto_iteration, true, true, false, num_viols);

  int buffer_moves_ = resizer_->buffer_move->numCommittedMoves();
  int size_up_moves_ = resizer_->size_up_move->numCommittedMoves();
  int size_down_moves_ = resizer_->size_down_move->numCommittedMoves();
  int swap_pins_moves_ = resizer_->swap_pins_move->numCommittedMoves();
  int clone_moves_ = resizer_->clone_move->numCommittedMoves();
  int split_load_moves_ = resizer_->split_load_move->numCommittedMoves();
  int unbuffer_moves_ = resizer_->unbuffer_move->numCommittedMoves();
  int ga_moves_ = resizer_->ga_move->numCommittedMoves();

  if (unbuffer_moves_ > 0) {
    repaired = true;
    logger_->info(RSZ, 59, "Removed {} buffers.", unbuffer_moves_);
  }
  if (buffer_moves_ > 0 || split_load_moves_ > 0) {
    repaired = true;
    if (split_load_moves_ == 0) {
      logger_->info(RSZ, 40, "Inserted {} buffers.", buffer_moves_);
    } else {
      logger_->info(RSZ,
                    45,
                    "Inserted {} buffers, {} to split loads.",
                    buffer_moves_ + split_load_moves_,
                    split_load_moves_);
    }
  }
  logger_->metric("design__instance__count__setup_buffer",
                  buffer_moves_ + split_load_moves_);
  if (size_up_moves_ + size_down_moves_ > 0) {
    repaired = true;
    logger_->info(RSZ,
                  51,
                  "Resized {} instances, {} sized up, {} sized down.",
                  size_up_moves_ + size_down_moves_,
                  size_up_moves_,
                  size_down_moves_);
  }
  if (swap_pins_moves_ > 0) {
    repaired = true;
    logger_->info(RSZ, 43, "Swapped pins on {} instances.", swap_pins_moves_);
  }
  if (clone_moves_ > 0) {
    repaired = true;
    logger_->info(RSZ, 49, "Cloned {} instances.", clone_moves_);
  }
  if (ga_moves_ > 0) {
    repaired = true;
    logger_->info(RSZ, 53, "Sized {} gates.", ga_moves_);
  }
  const Slack worst_slack = sta_->worstSlack(max_);
  if (fuzzyLess(worst_slack, setup_slack_margin)) {
    repaired = true;
    logger_->warn(RSZ, 62, "Unable to repair all setup violations.");
  }
  if (resizer_->overMaxArea()) {
    logger_->error(RSZ, 25, "max utilization reached.");
  }

  return repaired;
}

// For testing.
void RepairSetup::repairSetup(const Pin* end_pin)
{
  init();
  max_repairs_per_pass_ = 1;

  Vertex* vertex = graph_->pinLoadVertex(end_pin);
  const Slack slack = sta_->vertexSlack(vertex, max_);
  Path* path = sta_->vertexWorstSlackPath(vertex, max_);

  move_sequence.clear();
  move_sequence = {resizer_->unbuffer_move,
                   resizer_->size_down_move,
                   resizer_->size_up_move,
                   resizer_->swap_pins_move,
                   resizer_->buffer_move,
                   resizer_->clone_move,
                   resizer_->split_load_move};

  {
    IncrementalParasiticsGuard guard(resizer_);
    repairPath(path, slack, 0.0);
  }

  int unbuffer_moves_ = resizer_->unbuffer_move->numCommittedMoves();
  if (unbuffer_moves_ > 0) {
    logger_->info(RSZ, 61, "Removed {} buffers.", unbuffer_moves_);
  }
  int buffer_moves_ = resizer_->buffer_move->numCommittedMoves();
  int split_load_moves_ = resizer_->split_load_move->numMoves();
  if (buffer_moves_ + split_load_moves_ > 0) {
    logger_->info(
        RSZ, 30, "Inserted {} buffers.", buffer_moves_ + split_load_moves_);
  }
  int size_up_moves_ = resizer_->size_up_move->numMoves();
  int size_down_moves_ = resizer_->size_down_move->numMoves();
  if (size_up_moves_ + size_down_moves_ > 0) {
    logger_->info(RSZ,
                  38,
                  "Resized {} instances, {} sized up, {} sized down.",
                  size_up_moves_ + size_down_moves_,
                  size_up_moves_,
                  size_down_moves_);
  }
  int swap_pins_moves_ = resizer_->swap_pins_move->numMoves();
  if (swap_pins_moves_ > 0) {
    logger_->info(RSZ, 44, "Swapped pins on {} instances.", swap_pins_moves_);
  }
}

int RepairSetup::fanout(Vertex* vertex)
{
  int fanout = 0;
  VertexOutEdgeIterator edge_iter(vertex, graph_);
  while (edge_iter.hasNext()) {
    Edge* edge = edge_iter.next();
    // Disregard output->output timing arcs
    if (edge->isWire()) {
      fanout++;
    }
  }
  return fanout;
}

/* This is the main routine for repairing setup violations. We have
 - remove driver (step 1)
 - upsize driver (step 2)
 - rebuffer (step 3)
 - swap pin (step 4)
 - split loads (step 5)
 And they are always done in the same order. Not clear whether
 this order is the best way at all times. Also need to worry about
 actually using global routes...
 Things that can be added:
 - Intelligent rebuffering .... so if we added 2 buffers then maybe add
   two inverters instead.
 - pin swap (V0 is done)
 - Logic cloning
 - VT swap (already there via the normal resize code.... but we need to
   figure out how to deal with min implant rules to make it production
   ready)
 */
bool RepairSetup::repairPath(Path* path,
                             const Slack path_slack,
                             const float setup_slack_margin)
{
  // float initial_tns = sta_->totalNegativeSlack(max_);
  PathExpanded expanded(path, sta_);
  int changed = 0;

  if (expanded.size() > 1) {
    const int path_length = expanded.size();
    vector<pair<int, Delay>> load_delays;
    const int start_index = expanded.startIndex();
    const DcalcAnalysisPt* dcalc_ap = path->dcalcAnalysisPt(sta_);
    const int lib_ap = dcalc_ap->libertyIndex();
    // Find load delay for each gate in the path.
    for (int i = start_index; i < path_length; i++) {
      const Path* path = expanded.path(i);
      Vertex* path_vertex = path->vertex(sta_);
      const Pin* path_pin = path->pin(sta_);
      if (i > 0 && network_->isDriver(path_pin)
          && !network_->isTopLevelPort(path_pin)) {
        const TimingArc* prev_arc = path->prevArc(sta_);
        const TimingArc* corner_arc = prev_arc->cornerArc(lib_ap);
        Edge* prev_edge = path->prevEdge(sta_);
        const Delay load_delay
            = graph_->arcDelay(prev_edge, prev_arc, dcalc_ap->index())
              // Remove intrinsic delay to find load dependent delay.
              - corner_arc->intrinsicDelay();
        load_delays.emplace_back(i, load_delay);
        debugPrint(logger_,
                   RSZ,
                   "repair_setup",
                   3,
                   "{} load_delay = {} intrinsic_delay = {}",
                   path_vertex->name(network_),
                   delayAsString(load_delay, sta_, 3),
                   delayAsString(corner_arc->intrinsicDelay(), sta_, 3));
      }
    }

    sort(
        load_delays.begin(),
        load_delays.end(),
        [](pair<int, Delay> pair1, pair<int, Delay> pair2) {
          return pair1.second > pair2.second
                 || (pair1.second == pair2.second && pair1.first > pair2.first);
        });
    // Attack gates with largest load delays first.
    int repairs_per_pass = 1;
    if (max_viol_ - min_viol_ != 0.0) {
      repairs_per_pass
          += std::round((max_repairs_per_pass_ - 1) * (-path_slack - min_viol_)
                        / (max_viol_ - min_viol_));
    }
    if (fallback_) {
      repairs_per_pass = 1;
    }
    debugPrint(logger_,
               RSZ,
               "repair_setup",
               3,
               "Path slack: {}, repairs: {}, load_delays: {}",
               delayAsString(path_slack, sta_, 3),
               repairs_per_pass,
               load_delays.size());
    bool do_ga = true;
    for (const auto& [drvr_index, ignored] : load_delays) {
      if (changed >= repairs_per_pass) {
        break;
      }
      const Path* drvr_path = expanded.path(drvr_index);
      Vertex* drvr_vertex = drvr_path->vertex(sta_);
      const Pin* drvr_pin = drvr_vertex->pin();
      LibertyPort* drvr_port = network_->libertyPort(drvr_pin);
      LibertyCell* drvr_cell = drvr_port ? drvr_port->libertyCell() : nullptr;
      const int fanout = this->fanout(drvr_vertex);
      debugPrint(logger_,
                 RSZ,
                 "repair_setup",
                 3,
                 "{} {} fanout = {} drvr_index = {}",
                 network_->pathName(drvr_pin),
                 drvr_cell ? drvr_cell->name() : "none",
                 fanout,
                 drvr_index);

      for (BaseMove* move : move_sequence) {
        debugPrint(logger_,
                   RSZ,
                   "repair_setup",
                   1,
                   "Considering {} for {}",
                   move->name(),
                   network_->pathName(drvr_pin));
        if (!do_ga && ga_enabled_) {
          continue;
        }        
        if (strcmp(move->name(), "GAMove") == 0 && do_ga && ga_enabled_) {
          do_ga = false;
        }

        if (move->doMove(drvr_path,
                         drvr_index,
                         path_slack,
                         &expanded,
                         setup_slack_margin)) {
          if (move == resizer_->unbuffer_move) {
            // Only allow one unbuffer move per pass to
            // prevent the use-after-free error of multiple buffer removals.
            changed += repairs_per_pass;
          } else {
            changed++;
          }
          // Move on to the next gate
          break;
        }
        debugPrint(logger_,
                   RSZ,
                   "repair_setup",
                   2,
                   "Move {} failed for {}",
                   move->name(),
                   network_->pathName(drvr_pin));
      }
    }
  }
  // float curr_tns = sta_->totalNegativeSlack(max_);
  // std::cout << "curr_tns: " << curr_tns << " initial_tns: " << initial_tns << std::endl;
//   if (curr_tns < initial_tns) {
//     return true;
//   }
//   else {
//     return false;
//   }
  return changed > 0;
}

void RepairSetup::printProgress(const int iteration,
                                const bool force,
                                const bool end,
                                const bool last_gasp,
                                const int num_viols) const
{
  const bool start = iteration == 0;

  if (start && !end) {
    logger_->report(
        "   Iter   | Removed | Resized | Inserted | Cloned |  Pin  |GA Resized|"
        "   Area   |    WNS   |   TNS      |  Viol  | Worst");
    logger_->report(
        "          | Buffers |  Gates  | Buffers  |  Gates | Swaps |   Gates  |"
        "          |          |            | Endpts | Endpt");
    logger_->report(
        "-----------------------------------------------------------"
        "---------------------------------------------------");
  }

  if (iteration % print_interval_ == 0 || force || end) {
    Slack wns;
    Vertex* worst_vertex;
    // useful for CAD
    sta_->worstSlack(max_, wns, worst_vertex);
    const Slack tns = sta_->totalNegativeSlack(max_);

    std::string itr_field
        = fmt::format("{}{}", iteration, (last_gasp ? "*" : ""));
    if (end) {
      itr_field = "final";
    }

    const double design_area = resizer_->computeDesignArea();
    const double area_growth = design_area - initial_design_area_;

    // This actually prints both committed and pending moves, so the moves could
    // could go down if a pass is restrored by the journal.
    logger_->report(
        "{: >9s} | {: >7d} | {: >7d} | {: >8d} | {: >6d} | {: >5d} "
        "| {: >8d} | {: >+7.1f}% | {: >8s} | {: >10s} | {: >6d} | {}",
        itr_field,
        resizer_->unbuffer_move->numCommittedMoves(),
        resizer_->size_up_move->numCommittedMoves()
            + resizer_->size_down_move->numCommittedMoves(),
        resizer_->buffer_move->numCommittedMoves()
            + resizer_->split_load_move->numCommittedMoves(),
        resizer_->clone_move->numCommittedMoves(),
        resizer_->swap_pins_move->numCommittedMoves(),
        resizer_->ga_move->numCommittedMoves(),
        area_growth / initial_design_area_ * 1e2,
        delayAsString(wns, sta_, 3),
        delayAsString(tns, sta_, 1),
        max(0, num_viols),
        worst_vertex != nullptr ? worst_vertex->name(network_) : "");
  }

  if (end) {
    logger_->report(
        "-----------------------------------------------------------"
        "---------------------------------------------------");
  }
}

// Terminate progress if incremental fix rate within an opto interval falls
// below the threshold.   Bump up the threshold after each large opto
// interval.
bool RepairSetup::terminateProgress(const int iteration,
                                    const float initial_tns,
                                    float& prev_tns,
                                    float& fix_rate_threshold,
                                    // for debug only
                                    const int endpt_index,
                                    const int num_endpts)
{
  if (iteration % opto_large_interval_ == 0) {
    fix_rate_threshold *= 2.0;
  }
  if (iteration % opto_small_interval_ == 0) {
    float curr_tns = sta_->totalNegativeSlack(max_);
    float inc_fix_rate = (prev_tns - curr_tns) / initial_tns;
    prev_tns = curr_tns;
    if (iteration > 1000  // allow for some initial fixing for 1000 iterations
        && inc_fix_rate < fix_rate_threshold) {
      // clang-format off
      debugPrint(logger_, RSZ, "repair_setup", 1, "bailing out at iter {}"
                 " because incr fix rate {:0.2f}% is < {:0.2f}% [endpt {}/{}]",
                 iteration, inc_fix_rate*100, fix_rate_threshold*100,
                 endpt_index, num_endpts);
      // clang-format on
      return true;
    }
  }
  return false;
}

// Perform some last fixing based on sizing only.
// This is a greedy opto that does not degrade WNS or TNS.
// TODO: add VT swap
void RepairSetup::repairSetupLastGasp(const OptoParams& params, int& num_viols)
{
  // Sort remaining failing endpoints
  const VertexSet* endpoints = sta_->endpoints();
  vector<pair<Vertex*, Slack>> violating_ends;
  for (Vertex* end : *endpoints) {
    const Slack end_slack = sta_->vertexSlack(end, max_);
    if (end_slack < params.setup_slack_margin) {
      violating_ends.emplace_back(end, end_slack);
    }
  }
  std::stable_sort(violating_ends.begin(),
                   violating_ends.end(),
                   [](const auto& end_slack1, const auto& end_slack2) {
                     return end_slack1.second < end_slack2.second;
                   });
  num_viols = violating_ends.size();

  float curr_tns = sta_->totalNegativeSlack(max_);
  if (fuzzyGreaterEqual(curr_tns, 0)) {
    // clang-format off
    debugPrint(logger_, RSZ, "repair_setup", 1, "last gasp is bailing out "
               "because TNS is {:0.2f}", curr_tns);
    // clang-format on
    return;
  }

  // Don't do anything unless there was some progress from previous fixing
  if ((params.initial_tns - curr_tns) / params.initial_tns < 0.05) {
    // clang-format off
    debugPrint(logger_, RSZ, "repair_setup", 1, "last gasp is bailing out "
               "because TNS was reduced by < 5% from previous fixing");
    // clang-format on
    return;
  }

  int end_index = 0;
  int max_end_count = violating_ends.size();
  if (max_end_count == 0) {
    // clang-format off
    debugPrint(logger_, RSZ, "repair_setup", 1, "last gasp is bailing out "
               "because there are no violating endpoints");
    // clang-format on
    return;
  }
  // clang-format off
  debugPrint(logger_, RSZ, "repair_setup", 1, "{} violating endpoints remain",
             max_end_count);
  // clang-format on
  int opto_iteration = params.iteration;
  printProgress(opto_iteration, false, false, true, num_viols);

  float prev_tns = curr_tns;
  Slack curr_worst_slack = violating_ends[0].second;
  Slack prev_worst_slack = curr_worst_slack;
  bool prev_termination = false;
  bool two_cons_terminations = false;
  float fix_rate_threshold = inc_fix_rate_threshold_;

  for (const auto& end_original_slack : violating_ends) {
    fallback_ = false;
    Vertex* end = end_original_slack.first;
    Slack end_slack = sta_->vertexSlack(end, max_);
    Slack worst_slack;
    Vertex* worst_vertex;
    sta_->worstSlack(max_, worst_slack, worst_vertex);
    end_index++;
    if (end_index > max_end_count) {
      break;
    }
    int pass = 1;
    resizer_->journalBegin();
    while (pass <= max_last_gasp_passes_) {
      opto_iteration++;
      // if (terminateProgress(opto_iteration,
      //                       params.initial_tns,
      //                       prev_tns,
      //                       fix_rate_threshold,
      //                       end_index,
      //                       max_end_count)) {
      //   if (prev_termination) {
      //     // Abort entire fixing if no progress for 200 iterations
      //     two_cons_terminations = true;
      //   } else {
      //     prev_termination = true;
      //   }
      //   resizer_->journalEnd();
      //   break;
      // }
      if (opto_iteration % opto_small_interval_ == 0) {
        prev_termination = false;
      }
      if (params.verbose || opto_iteration == 1) {
        printProgress(opto_iteration, false, false, true, num_viols);
      }
      if (end_slack > params.setup_slack_margin) {
        --num_viols;
        resizer_->journalEnd();
        break;
      }
      Path* end_path = sta_->vertexWorstSlackPath(end, max_);

      const bool changed
          = repairPath(end_path, end_slack, params.setup_slack_margin);

      if (!changed) {
        if (pass != 1) {
          resizer_->journalRestore();
        } else {
          resizer_->journalEnd();
        }
        break;
      }
      resizer_->updateParasitics();
      sta_->findRequireds();
      end_slack = sta_->vertexSlack(end, max_);
      sta_->worstSlack(max_, curr_worst_slack, worst_vertex);
      curr_tns = sta_->totalNegativeSlack(max_);

      // Accept only moves that improve both WNS and TNS
      if (fuzzyGreaterEqual(curr_tns, prev_tns)) {
        // clang-format off
        debugPrint(logger_, RSZ, "repair_setup", 1, "sizing move accepted for "
                   "endpoint {} pass {} because WNS improved to {:0.3f} and "
                   "TNS improved to {:0.3f}",
                   end_index, pass, curr_worst_slack, curr_tns);
        // clang-format on
        prev_worst_slack = curr_worst_slack;
        prev_tns = curr_tns;
        if (end_slack > params.setup_slack_margin) {
          --num_viols;
        }
        resizer_->journalEnd();
        resizer_->journalBegin();
      } else {
        fallback_ = true;
        resizer_->journalRestore();
        break;
      }

      if (resizer_->overMaxArea()) {
        resizer_->journalEnd();
        break;
      }
      if (end_index == 1) {
        end = worst_vertex;
      }
      pass++;
    }  // while pass <= max_last_gasp_passes_
    if (params.verbose || opto_iteration == 1) {
      printProgress(opto_iteration, true, false, true, num_viols);
    }
    if (two_cons_terminations) {
      // clang-format off
      debugPrint(logger_, RSZ, "repair_setup", 1, "bailing out of last gasp fixing"
                 "due to no TNS progress for two opto cycles");
      // clang-format on
      break;
    }
  }  // for each violating endpoint
}

// Genetic Algorithm implementation for gate sizing optimization
bool RepairSetup::gateSizingWithGa(const rsz::GaParams& ga_params,
                                     const float setup_slack_margin,
                                     const bool verbose)
{
  init();
  logger_->info(RSZ, 206, "Starting Genetic Algorithm for gate sizing optimization");
  
  // Initialize random generator
  std::random_device rd;
  ga_random_generator_.seed(rd());
  
  // Initialize gates for GA optimization
  initializeGaGates();
  
  if (ga_gates_.empty()) {
    logger_->info(RSZ, 207, "No gates available for GA optimization");
    return false;
  }
  
  logger_->info(RSZ, 208, "GA optimizing {} gates with {} generations, population size {}",
                ga_gates_.size(), ga_params.max_generations, ga_params.population_size);
  
  // Store original design state
  restoreOriginalSizes();  // This stores the current state as original
  
  // Initialize population
  std::vector<GaChromosome> population(ga_params.population_size);
  initializeGaPopulation(population, ga_params);
  GaChromosome chromosome;
  for (size_t i = 0; i < ga_gates_.size(); i++) {
    chromosome.genes.resize(ga_gates_.size());
    chromosome.genes[i] = ga_gates_[i].current_size_index;
  }
  population.push_back(chromosome);
  
  GaChromosome best_overall;
  best_overall.fitness = -std::numeric_limits<float>::max();
  
  IncrementalParasiticsGuard guard(resizer_);
  
  // Evolution loop
  for (int generation = 0; generation < ga_params.max_generations; generation++) {
    // Evaluate fitness for all chromosomes
    for (auto& chromosome : population) {
      chromosome.fitness = evaluateGaFitness(chromosome);
    }
    
    // Sort population by fitness (best first)
    std::sort(population.begin(), population.end());
    
    // Update best overall
    if (population[0].fitness > best_overall.fitness) {
      best_overall = population[0];
    }
    
    // Print progress
    if (verbose || generation % 10 == 0) {
      printGaProgress(generation, population[0], verbose);
    }
    
    // Check termination conditions
    if (population[0].wns >= setup_slack_margin) {
      logger_->info(RSZ, 209, "GA converged at generation {} with WNS = {}",
                    generation, population[0].wns);
      break;
    }
    
    // Create new generation
    std::vector<GaChromosome> new_population;
    
    // Elitism: keep best chromosomes
    int elite_count = static_cast<int>(ga_params.population_size * ga_params.elitism_rate);
    for (int i = 0; i < elite_count; i++) {
      new_population.push_back(population[i]);
    }
    
    // Generate offspring through crossover and mutation
    while (new_population.size() < ga_params.population_size) {
      GaChromosome parent1 = gaSelection(population);
      GaChromosome parent2 = gaSelection(population);
      
      GaChromosome offspring = gaCrossover(parent1, parent2, ga_params.crossover_rate);
      gaMutation(offspring, ga_params.mutation_rate);
      
      new_population.push_back(offspring);
    }
    
    population = std::move(new_population);
  }
  
  // Apply best solution
  applyGaChromosome(best_overall);
  resizer_->updateParasitics();
  sta_->findRequireds();
  
  logger_->info(RSZ, 210, "GA optimization completed. Best fitness: {}, WNS: {}, TNS: {}",
                best_overall.fitness, best_overall.wns, best_overall.tns);
  
  return best_overall.wns >= setup_slack_margin;
}

void RepairSetup::initializeGaGates()
{
  ga_gates_.clear();
  original_cells_.clear();
  
  // Collect all resizable gates in the design
  Instance* top_inst = network_->topInstance();
  // std::cout << "instance count: " << network_->instanceCount() << std::endl;
  
  // Method 1: Use leafInstanceIterator for all leaf instances
  sta::InstanceChildIterator *child_iter = network_->childIterator(top_inst);
  while (child_iter->hasNext()) {
    Instance *inst = child_iter->next();
    // std::cout << "child instance found: " << network_->name(inst) << std::endl;
    // std::cout << "liberty cell found: " << network_->cellName(inst) << std::endl;
    
    LibertyCell* lib_cell = network_->libertyCell(inst);
    if (!lib_cell) {
      continue;
    }
    
    // Find equivalent cells (different sizes of the same functionality)
    std::vector<LibertyCell*> equiv_cells;
    
    // TODO: Add logic to find equivalent cells with different sizes
    // This would typically involve checking the cell library for cells with
    // the same function but different sizes/drives
    int current_size_index = 0;
    LibertyCellSeq swappable_cells = resizer_->getSwappableCells(lib_cell);
    sort(swappable_cells,
         [=](const LibertyCell* cell1, const LibertyCell* cell2) {
           // Calculate metrics for cell1
           float input_cap1 = 0.0f, output_drive_res1 = FLT_MAX;
           float output_intrinsic_delay1 = 0.0f;
           sta::LibertyCellPortIterator port_iter1(cell1);
           while (port_iter1.hasNext()) {
             LibertyPort* port = port_iter1.next();
             if (port->direction()->isInput()) {
               input_cap1 += port->capacitance();
             } else if (port->direction()->isOutput()) {
               output_drive_res1 = std::min(output_drive_res1, port->driveResistance());
               output_intrinsic_delay1 += port->intrinsicDelay(this);
             }
           }
           
           // Calculate metrics for cell2
           float input_cap2 = 0.0f, output_drive_res2 = FLT_MAX;
           float output_intrinsic_delay2 = 0.0f;
           sta::LibertyCellPortIterator port_iter2(cell2);
           while (port_iter2.hasNext()) {
             LibertyPort* port = port_iter2.next();
             if (port->direction()->isInput()) {
               input_cap2 += port->capacitance();
             } else if (port->direction()->isOutput()) {
               output_drive_res2 = std::min(output_drive_res2, port->driveResistance());
               output_intrinsic_delay2 += port->intrinsicDelay(this);
             }
           }
           
            // Sort by: 1) Stronger drive (smaller resistance), 2) Faster (smaller delay), 3) Smaller input cap
            return std::tie(output_drive_res1, output_intrinsic_delay1, input_cap1)
                   < std::tie(output_drive_res2, output_intrinsic_delay2, input_cap2);
         });
    
    for (LibertyCell* cell : swappable_cells) {
      equiv_cells.push_back(cell);
      if (cell->name() == lib_cell->name()) {
        current_size_index = equiv_cells.size() - 1;
      }
      // std::cout << "swappable cell found: " << cell->name() << std::endl;
    }
    
    if (equiv_cells.size() > 1) {  // Only include if multiple sizes available
      GaGateInfo gate_info(inst);
      gate_info.available_cells = equiv_cells;
      gate_info.current_size_index = current_size_index;  // Assume current cell is at index 0
      ga_gates_.push_back(gate_info);
      original_cells_.push_back(lib_cell);
    }
  }
  delete child_iter;
  /*
  Alternative Methods to Iterate Instances:
  
  // Method 2: Use leafInstances() to get all leaf instances as a sequence
  InstanceSeq leaf_instances = network_->leafInstances();
  for (Instance* inst : leaf_instances) {
    // Process each instance...
  }
  
  // Method 3: Manually iterate children of top instance (hierarchical)
  InstanceChildIterator* child_iter = network_->childIterator(network_->topInstance());
  while (child_iter->hasNext()) {
    Instance* inst = child_iter->next();
    // Process each child instance...
    // Note: This only gets direct children, not all descendants
  }
  delete child_iter;
  
  // Method 4: Recursive function to traverse all instances in hierarchy
  void traverseInstances(Instance* parent) {
    InstanceChildIterator* child_iter = network_->childIterator(parent);
    while (child_iter->hasNext()) {
      Instance* inst = child_iter->next();
      // Process instance...
      if (network_->isHierarchical(inst)) {
        traverseInstances(inst);  // Recursive call for hierarchical instances
      }
    }
    delete child_iter;
  }
  */
}

void RepairSetup::initializeGaPopulation(std::vector<GaChromosome>& population,
                                          const GaParams& params)
{
  std::uniform_int_distribution<int> gene_dist;
  
  for (auto& chromosome : population) {
    chromosome.genes.resize(ga_gates_.size());
    
    for (size_t i = 0; i < ga_gates_.size(); i++) {
      // Random size selection for each gate
      int max_size_index = ga_gates_[i].available_cells.size() - 1;
      std::uniform_int_distribution<int> size_dist(0, max_size_index);
      chromosome.genes[i] = size_dist(ga_random_generator_);
    }
  }
}

float RepairSetup::evaluateGaFitness(const GaChromosome& chromosome)
{
  // Apply chromosome to design
  for (size_t i = 0; i < ga_gates_.size(); i++) {
    Instance* inst = ga_gates_[i].instance;
    LibertyCell* new_cell = ga_gates_[i].available_cells[chromosome.genes[i]];
    
    // Replace cell (this is a simplified version - actual implementation
    // would need to handle pin mapping, etc.)
    resizer_->replaceCell(inst, new_cell, true);
  }
  
  // Update parasitic and timing
  resizer_->updateParasitics();
  sta_->findRequireds();
  
  // Calculate metrics
  Slack wns;
  Vertex* worst_vertex;
  sta_->worstSlack(max_, wns, worst_vertex);
  Slack tns = sta_->totalNegativeSlack(max_);
  
  // Store metrics in chromosome (const_cast for this evaluation)
  const_cast<GaChromosome&>(chromosome).wns = wns;
  const_cast<GaChromosome&>(chromosome).tns = tns;

  // logger_->info(RSZ, 211, "GA evaluation: WNS: {}, TNS: {}", wns, tns);
  
  // Calculate fitness function
  // Higher fitness is better
  // Prioritize timing improvement, then minimize area
  float fitness = 0.0f;
  
  // Timing violations exist - improve timing
  fitness = wns * 100.0f + (tns * 10.0f);
  
  // logger_->info(RSZ, 212, "GA fitness: {}", fitness);
  
  return fitness;
}

void RepairSetup::applyGaChromosome(const GaChromosome& chromosome)
{
  for (size_t i = 0; i < ga_gates_.size(); i++) {
    Instance* inst = ga_gates_[i].instance;
    LibertyCell* new_cell = ga_gates_[i].available_cells[chromosome.genes[i]];
    resizer_->replaceCell(inst, new_cell, true);
  }
}

void RepairSetup::restoreOriginalSizes()
{
  for (size_t i = 0; i < ga_gates_.size(); i++) {
    Instance* inst = ga_gates_[i].instance;
    LibertyCell* original_cell = original_cells_[i];
    resizer_->replaceCell(inst, original_cell, true);
  }
}

GaChromosome RepairSetup::gaSelection(const std::vector<GaChromosome>& population)
{
  // Tournament selection
  const int tournament_size = 3;
  std::uniform_int_distribution<int> pop_dist(0, population.size() - 1);
  
  GaChromosome best = population[pop_dist(ga_random_generator_)];
  
  for (int i = 1; i < tournament_size; i++) {
    GaChromosome candidate = population[pop_dist(ga_random_generator_)];
    if (candidate.fitness > best.fitness) {
      best = candidate;
    }
  }
  
  return best;
}

GaChromosome RepairSetup::gaCrossover(const GaChromosome& parent1,
                                       const GaChromosome& parent2,
                                       const float crossover_rate)
{
  GaChromosome offspring;
  offspring.genes.resize(parent1.genes.size());
  
  std::uniform_real_distribution<float> real_dist(0.0f, 1.0f);
  
  if (real_dist(ga_random_generator_) < crossover_rate) {
    // Single-point crossover
    // std::uniform_int_distribution<int> point_dist(1, parent1.genes.size() - 1);
    // int crossover_point = point_dist(ga_random_generator_);
    
    // for (size_t i = 0; i < offspring.genes.size(); i++) {
    //   if (i < crossover_point) {
    //     offspring.genes[i] = parent1.genes[i];
    //   } else {
    //     offspring.genes[i] = parent2.genes[i];
    //   }
    // }

    // probability of selecting crossover is based on fitness
    float selection_rate = parent2.fitness / (parent1.fitness + parent2.fitness);
    for (size_t i = 0; i < offspring.genes.size(); i++) {
      if (real_dist(ga_random_generator_) < selection_rate) {
        offspring.genes[i] = parent1.genes[i];
      } else {
        offspring.genes[i] = parent2.genes[i];
      }
    }
  } else {
    // No crossover - copy parent1
    offspring.genes = parent1.genes;
  }
  
  return offspring;
}

void RepairSetup::gaMutation(GaChromosome& chromosome, const float mutation_rate)
{
  std::uniform_real_distribution<float> real_dist(0.0f, 1.0f);
  
  for (size_t i = 0; i < chromosome.genes.size(); i++) {
    if (real_dist(ga_random_generator_) < mutation_rate) {
      // Mutate this gene
      int max_size_index = ga_gates_[i].available_cells.size() - 1;
      std::uniform_int_distribution<int> size_dist(0, max_size_index);
      chromosome.genes[i] = size_dist(ga_random_generator_);
    }
  }
}

void RepairSetup::printGaProgress(const int generation,
                                  const GaChromosome& best_chromosome,
                                  const bool verbose) const
{
  if (generation == 0) {
    logger_->report("Gen |   Fitness   |    WNS     |    TNS     ");
    logger_->report("----+-------------+------------+------------");
  }
  
  logger_->report("{:3d} | {: >10s} | {: >8s} | {: >10s}",
                  generation,
                  delayAsString(best_chromosome.fitness, sta_, 3),
                  delayAsString(best_chromosome.wns, sta_, 3),
                  delayAsString(best_chromosome.tns, sta_, 1));
}
}  // namespace rsz
