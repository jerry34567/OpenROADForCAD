// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2025-2025, The OpenROAD Authors

#include "GAMove.hh"

#include <cmath>
#include <algorithm>

#include "BaseMove.hh"

namespace rsz {

using std::string;

using utl::RSZ;

using sta::ArcDelay;
using sta::DcalcAnalysisPt;
using sta::Instance;
using sta::InstancePinIterator;
using sta::LibertyCell;
using sta::LibertyPort;
using sta::LoadPinIndexMap;
using sta::NetConnectedPinIterator;
using sta::Path;
using sta::PathExpanded;
using sta::Pin;
using sta::Slack;
using sta::Slew;

// Calculate total delay of all stages in the path
float GAMove::calculateTotalPathDelay(PathExpanded* expanded, 
                                          int start_index, 
                                          int end_index)
{
    float total_delay = 0.0f;
    
    // Iterate through each pin in the path
    for (int i = start_index; i < end_index; i += 2) {
        const Path* drvr_path = expanded->path(i);
        Pin* drvr_pin = drvr_path->pin(sta_);
        const DcalcAnalysisPt* dcalc_ap = drvr_path->dcalcAnalysisPt(sta_);
        const int lib_ap = dcalc_ap->libertyIndex();
        
        // Get driver port
        LibertyPort* drvr_port = network_->libertyPort(drvr_pin);
        if (!drvr_port) {
            continue; // Skip if no liberty port
        }
        
        // Get gate instance to calculate intrinsic delay
        Instance* gate_inst = network_->instance(drvr_pin);
        LibertyCell* gate_cell = network_->libertyCell(gate_inst);
        
        // Calculate intrinsic delay of the gate
        float intrinsic_delay = 0.0f;
        if (gate_cell) {
            // Get intrinsic delay from the output port
            intrinsic_delay = drvr_port->intrinsicDelay(sta_);
        }
        
        // Get drive resistance of output port
        const float drive_resistance = drvr_port->driveResistance();
        
        // Calculate total fanout capacitance
        float total_fanout_cap = 0.0f;
        
        // Get the net connected to this driver pin
        Net* drvr_net = network_->net(drvr_pin);
        if (drvr_net) {
            // Iterate through all pins connected to this net
            NetConnectedPinIterator* pin_iter = network_->connectedPinIterator(drvr_net);
            while (pin_iter->hasNext()) {
                const Pin* connected_pin = pin_iter->next();
                
                // Skip the driver pin itself
                if (connected_pin == drvr_pin) {
                    continue;
                }
                
                // Add capacitance of load pins
                LibertyPort* connected_port = network_->libertyPort(connected_pin);
                if (connected_port && connected_port->direction()->isInput()) {
                    total_fanout_cap += connected_port->cornerPort(lib_ap)->capacitance();
                }
            }
            delete pin_iter;
        }
        
        // Calculate gate contribution: drive_resistance * total_fanout_cap + intrinsic_delay
        const float gate_contribution = drive_resistance * total_fanout_cap + intrinsic_delay;
        total_delay += gate_contribution;
        
        debugPrint(logger_,
                   RSZ,
                   "path_delay",
                   3,
                   "Stage {}: drive_resistance={}, total_fanout_cap={}, intrinsic_delay={}, gate_contribution={}",
                   (i/2),
                   drive_resistance,
                   total_fanout_cap,
                   delayAsString(intrinsic_delay, sta_, 3),
                   delayAsString(gate_contribution, sta_, 3));
    }
    
    debugPrint(logger_,
               RSZ,
               "path_delay",
               2,
               "Total path delay: {}",
               delayAsString(total_delay, sta_, 3));
    
    return total_delay;
}

bool GAMove::doMove(const Path* drvr_path,
                        int drvr_index,
                        Slack drvr_slack,
                        PathExpanded* expanded,
                        float setup_slack_margin)
{
  path_length_ = expanded->size();
  start_index_ = expanded->startIndex();
  expanded_ = expanded;
  std::vector<Instance*> insts;
  // Find each gate in the path.
  for (int i = start_index_; i < path_length_; i++) {
    const Path* path = expanded->path(i);
    const Pin* path_pin = path->pin(sta_);
    Instance* inst = network_->instance(path_pin);
    insts.push_back(inst);
  }

  logger_->info(RSZ, 213, "GAMove::doMove size {}", insts.size());

  // Initialize random generator
  std::random_device rd;
  ga_random_generator_.seed(rd());

  logger_->info(RSZ, 214, "Starting Genetic Algorithm for gate sizing optimization");
  
  
  
  // Initialize gates for GA optimization
  initializeGaGates(insts);
  
  if (ga_gates_.empty()) {
    logger_->info(RSZ, 215, "No gates available for GA optimization");
    return false;
  }

  GaParams ga_params;
  ga_params.max_generations = 30;
  ga_params.population_size = 50;
  ga_params.elitism_rate = 0.2;
  ga_params.crossover_rate = 0.8;
  ga_params.mutation_rate = 0.1;
  
  logger_->info(RSZ, 216, "GA optimizing {} gates with {} generations, population size {}",
                ga_gates_.size(), ga_params.max_generations, ga_params.population_size);
  
  // Store original design state
  // restoreOriginalSizes();  // This stores the current state as original
  
  // Initialize population
  std::vector<GaChromosome> population(ga_params.population_size);
  initializeGaPopulation(population, ga_params);
  GaChromosome chromosome;
  for (size_t i = 0; i < ga_gates_.size(); i++) {
    chromosome.genes.resize(ga_gates_.size());
    chromosome.genes[i] = ga_gates_[i].current_size_index;
  }
  population.push_back(chromosome);
  logger_->info(RSZ, 218, "Initial fitness: {}", evaluateGaFitness(chromosome));
  
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
    if (generation % 10 == 0) {
      printGaProgress(generation, population[0], true);
    }
    // Check termination conditions
    // if (population[0].wns >= setup_slack_margin) {
    //   logger_->info(RSZ, 217, "GA converged at generation {} with WNS = {}",
    //                 generation, population[0].wns);
    //   break;
    // }
    
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

  bool changed = false;
  for (size_t i = 0; i < ga_gates_.size(); i++) {
    Instance* inst = ga_gates_[i].instance;
    LibertyCell* new_cell = ga_gates_[i].available_cells[best_overall.genes[i]];
    if (new_cell != original_cells_[i]) {
      changed = true;
      addMove(inst, 1);
    }
  }

  return changed;
}

void GAMove::initializeGaGates(std::vector<Instance*> insts)
{
  ga_gates_.clear();
  original_cells_.clear();
  
  // Collect all resizable gates in the design
  for (Instance* inst : insts) {
    
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
                   > std::tie(output_drive_res2, output_intrinsic_delay2, input_cap2);
         });
    
    for (LibertyCell* cell : swappable_cells) {
      equiv_cells.push_back(cell);
      if (cell->name() == lib_cell->name()) {
        current_size_index = equiv_cells.size() - 1;
      }
    }
    
    if (equiv_cells.size() > 1) {  // Only include if multiple sizes available
      GaGateInfo gate_info(inst);
      gate_info.available_cells = equiv_cells;
      gate_info.current_size_index = current_size_index;  // Assume current cell is at index 0
      ga_gates_.push_back(gate_info);
      original_cells_.push_back(lib_cell);
    }
  }
}

void GAMove::initializeGaPopulation(std::vector<GAMove::GaChromosome>& population,
                                          const GaParams& params)
{
  std::uniform_int_distribution<int> gene_dist;
  
  for (auto& chromosome : population) {
    chromosome.genes.resize(ga_gates_.size());
    
    for (size_t i = 0; i < ga_gates_.size(); i++) {
      // Random size selection for each gate
      int max_size_index = ga_gates_[i].available_cells.size() - 1;
      std::uniform_int_distribution<int> size_dist(std::max(0, ga_gates_[i].current_size_index - 1), max_size_index);
      chromosome.genes[i] = size_dist(ga_random_generator_);
    }
  }
}

float GAMove::evaluateGaFitness(const GAMove::GaChromosome& chromosome)
{
  // Apply chromosome to design
  for (size_t i = 0; i < ga_gates_.size(); i++) {
    Instance* inst = ga_gates_[i].instance;
    LibertyCell* new_cell = ga_gates_[i].available_cells[chromosome.genes[i]];
    
    // Replace cell (this is a simplified version - actual implementation
    // would need to handle pin mapping, etc.)
    resizer_->replaceCell(inst, new_cell, true);
  }

  // float fitness = calculateTotalPathDelay(expanded_, start_index_, path_length_);

  // // Update parasitic and timing
  resizer_->updateParasitics();
  sta_->findRequireds();
  const MinMax* max_ = MinMax::max();
  // Vertex* end = expanded_->path(path_length_ - 1)->vertex(sta_);
  // Slack end_slack = sta_->vertexSlack(end, max_);
  Slack tns = sta_->totalNegativeSlack(max_);
  float fitness = tns;
  
  // // Calculate metrics
  // Slack wns;
  // Vertex* worst_vertex;
  // const MinMax* max_ = MinMax::max();
  // sta_->worstSlack(max_, wns, worst_vertex);
  // Slack tns = sta_->totalNegativeSlack(max_);
  
  // Timing violations exist - improve timing
  // fitness = wns * 100.0f + (tns * 10.0f);
  
  return fitness;
}

void GAMove::applyGaChromosome(const GAMove::GaChromosome& chromosome)
{
  for (size_t i = 0; i < ga_gates_.size(); i++) {
    Instance* inst = ga_gates_[i].instance;
    LibertyCell* new_cell = ga_gates_[i].available_cells[chromosome.genes[i]];
    resizer_->replaceCell(inst, new_cell, true);
  }
}

void GAMove::restoreOriginalSizes()
{
  for (size_t i = 0; i < ga_gates_.size(); i++) {
    Instance* inst = ga_gates_[i].instance;
    LibertyCell* original_cell = original_cells_[i];
    resizer_->replaceCell(inst, original_cell, true);
  }
}

GAMove::GaChromosome GAMove::gaSelection(const std::vector<GaChromosome>& population)
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

GAMove::GaChromosome GAMove::gaCrossover(const GaChromosome& parent1,
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

void GAMove::gaMutation(GaChromosome& chromosome, const float mutation_rate)
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

void GAMove::printGaProgress(const int generation,
                                  const GaChromosome& best_chromosome,
                                  const bool verbose) const
{
  if (generation == 0) {
    logger_->report("Gen |   Fitness    ");
    logger_->report("----+-------------");
  }
  
  logger_->report("{:3d} | {: >10s}",
                  generation,
                  delayAsString(best_chromosome.fitness, sta_, 3));
}

// namespace rsz
}  // namespace rsz
