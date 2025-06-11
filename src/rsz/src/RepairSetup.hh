// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2022-2025, The OpenROAD Authors

#pragma once
#include <boost/functional/hash.hpp>
#include <queue>
#include <random>
#include <unordered_set>
#include <vector>

#include "db_sta/dbNetwork.hh"
#include "db_sta/dbSta.hh"
#include "rsz/Resizer.hh"
#include "sta/FuncExpr.hh"
#include "sta/MinMax.hh"
#include "sta/StaState.hh"
#include "utl/Logger.h"

namespace sta {
class PathExpanded;
}

namespace rsz {

class Resizer;
class RemoveBuffer;
class BaseMove;

using odb::Point;
using utl::Logger;

using sta::Corner;
using sta::dbNetwork;
using sta::dbSta;
using sta::DcalcAnalysisPt;
using sta::Delay;
using sta::Instance;
using sta::LibertyCell;
using sta::LibertyPort;
using sta::MinMax;
using sta::Net;
using sta::Path;
using sta::PathExpanded;
using sta::Pin;
using sta::RiseFall;
using sta::RiseFallBoth;
using sta::Slack;
using sta::Slew;
using sta::StaState;
using sta::TimingArc;
using sta::Vertex;

struct OptoParams
{
  int iteration;
  float initial_tns;
  const float setup_slack_margin;
  const bool verbose;

  OptoParams(const float margin, const bool verbose)
      : setup_slack_margin(margin), verbose(verbose)
  {
    iteration = 0;
    initial_tns = 0.0;
  }
};

struct GaChromosome
{
  std::vector<int> genes;  // Each gene represents the cell size index for a gate
  float fitness;
  float wns;  // Worst Negative Slack
  float tns;  // Total Negative Slack
  
  GaChromosome() : fitness(0.0f), wns(0.0f), tns(0.0f) {}
  
  bool operator<(const GaChromosome& other) const {
    return fitness > other.fitness;  // Higher fitness is better
  }
};

struct GaGateInfo
{
  Instance* instance;
  std::vector<LibertyCell*> available_cells;  // Available cell sizes
  int current_size_index;
  
  GaGateInfo(Instance* inst) : instance(inst), current_size_index(0) {}
};

class RepairSetup : public sta::dbStaState
{
 public:
  RepairSetup(Resizer* resizer);
  bool repairSetup(float setup_slack_margin,
                   // Percent of violating ends to repair to
                   // reduce tns (0.0-1.0).
                   double repair_tns_end_percent,
                   int max_passes,
                   int max_repairs_per_pass,
                   bool verbose,
                   const std::vector<MoveType>& sequence,
                   bool skip_pin_swap,
                   bool skip_gate_cloning,
                   bool skip_gate_sizing,
                   bool skip_size_down,
                   bool skip_buffering,
                   bool skip_split_load,
                   bool skip_buffer_removal,
                   bool skip_last_gasp,
                   bool ga_enabled,
                   bool shuffle_enabled);

  // Genetic Algorithm for gate sizing optimization
  bool gateSizingWithGa(const GaParams& ga_params,
                        const float setup_slack_margin,
                        const bool verbose);
  
  // For testing.
  void repairSetup(const Pin* end_pin);
  // For testing.
  void reportSwappablePins();
  // Rebuffer one net (for testing).
  // resizerPreamble() required.

 private:
  void init();
  bool repairPath(Path* path, Slack path_slack, float setup_slack_margin);
  int fanout(Vertex* vertex);
  bool hasTopLevelOutputPort(Net* net);

  void printProgress(int iteration,
                     bool force,
                     bool end,
                     bool last_gasp,
                     int num_viols) const;
  bool terminateProgress(int iteration,
                         float initial_tns,
                         float& prev_tns,
                         float& fix_rate_threshold,
                         int endpt_index,
                         int num_endpts);
  void repairSetupLastGasp(const OptoParams& params, int& num_viols);

  // Genetic Algorithm helper functions
  void initializeGaGates();
  void initializeGaPopulation(std::vector<GaChromosome>& population,
                              const GaParams& params);
  float evaluateGaFitness(const GaChromosome& chromosome);
  void applyGaChromosome(const GaChromosome& chromosome);
  void restoreOriginalSizes();
  GaChromosome gaSelection(const std::vector<GaChromosome>& population);
  GaChromosome gaCrossover(const GaChromosome& parent1, 
                           const GaChromosome& parent2,
                           float crossover_rate);
  void gaMutation(GaChromosome& chromosome, float mutation_rate);
  void printGaProgress(int generation, 
                       const GaChromosome& best_chromosome,
                       bool verbose) const;

  Logger* logger_ = nullptr;
  dbNetwork* db_network_ = nullptr;
  Resizer* resizer_;

  bool fallback_ = false;
  float min_viol_ = 0.0;
  float max_viol_ = 0.0;
  int max_repairs_per_pass_ = 1;
  int removed_buffer_count_ = 0;
  double initial_design_area_ = 0;

  std::vector<BaseMove*> move_sequence;

  bool ga_enabled_ = false;
  bool shuffle_enabled_ = false;

  // Genetic Algorithm data members
  std::vector<GaGateInfo> ga_gates_;
  std::vector<LibertyCell*> original_cells_;  // Store original cell assignments
  std::mt19937 ga_random_generator_;

  std::queue<Delay> shuffle_tns_windows_;
  float shuffle_tns_slope_ = 0.05;
  float shuffle_percent_ = 0.01;
  float best_tns_ = 0.0;
  int shuffle_times_ = 10;

  const MinMax* min_ = MinMax::min();
  const MinMax* max_ = MinMax::max();

  sta::UnorderedMap<LibertyPort*, sta::LibertyPortSet> equiv_pin_map_;

  static constexpr int decreasing_slack_max_passes_ = 50;
  static constexpr int print_interval_ = 10;
  static constexpr int opto_small_interval_ = 100;
  static constexpr int opto_large_interval_ = 1000;
  static constexpr float inc_fix_rate_threshold_
      = 0.0001;  // default fix rate threshold = 0.01%
  static constexpr int max_last_gasp_passes_ = 10;
};

}  // namespace rsz
