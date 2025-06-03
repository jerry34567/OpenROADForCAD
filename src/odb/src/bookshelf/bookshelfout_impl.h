// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2019-2025, The OpenROAD Authors

#pragma once

#include <list>
#include <map>
#include <string>

#include "odb/db.h"
#include "odb/dbMap.h"
#include "odb/bookshelf.h"
#include "odb/odb.h"
namespace utl {
class Logger;
}

namespace odb {

class dbBlock;
class dbBTerm;
class dbInst;
class dbTechNonDefaultRule;
class dbTechLayerRule;

class bookshelfout_impl
{
  FILE* _out;
  utl::Logger* _logger;
  double _dist_factor;

  int defdist(int value) { return (int) (((double) value) * _dist_factor); }
  int defdist(uint value) { return (uint) (((double) value) * _dist_factor); }

public:
  bookshelfout_impl(utl::Logger* logger) {
    _out = nullptr;
    _logger = logger;
    _dist_factor = 0;
  };
  ~bookshelfout_impl() = default;
  void writeInst(dbInst* inst);
  void writeBTerm(dbBTerm* bterm);
  void writeBPin(dbBPin* bpin, int cnt);
  bool writeBlock(dbBlock* block, const char* filename);
};

}  // namespace odb
