// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2025-2025, The OpenROAD Authors

#include <cmath>
#include <random>
#include <vector>
#include <map>

#include "BaseMove.hh"

namespace rsz {

// Forward declaration
struct GaParams;


class GAMove : public BaseMove
{
 public:
  using BaseMove::BaseMove;

  struct GaGateInfo
  {
    Instance* instance;
    std::vector<LibertyCell*> available_cells;  // Available cell sizes
    int current_size_index;
    
    GaGateInfo(Instance* inst) : instance(inst), current_size_index(0) {}
  };

  struct GaChromosome
  {
    std::vector<int> genes;  // Each gene represents the cell size index for a gate
    float fitness;
    
    GaChromosome() : fitness(0.0f) {}
    
    bool operator<(const GaChromosome& other) const {
      return fitness > other.fitness;  // Higher fitness is better
    }
  };

  bool doMove(const Path* drvr_path,
              int drvr_index,
              Slack drvr_slack,
              PathExpanded* expanded,
              float setup_slack_margin) override;

  const char* name() override { return "GAMove"; }

  // Genetic Algorithm helper functions
  void initializeGaGates(std::vector<Instance*> insts);
  void initializeGaPopulation(std::vector<GaChromosome>& population,
                              const GaParams& params);
  float evaluateGaFitness(const GaChromosome& chromosome);
  float evaluateGaFitnessAccurate(const GaChromosome& chromosome);
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
  
  // Genetic Algorithm data members
  std::vector<GaGateInfo> ga_gates_;
  std::vector<LibertyCell*> original_cells_;  // Store original cell assignments
  std::mt19937 ga_random_generator_;
  
  // Fitness evaluation cache to avoid redundant calculations
  std::map<std::vector<int>, float> fitness_cache_;
  int start_index_;
  int path_length_;
  PathExpanded* expanded_;
  
  // Calculate total delay of all stages in the path
  float calculateTotalPathDelay(PathExpanded* expanded, 
                                int start_index, 
                                int end_index);
};

}  // namespace rsz
