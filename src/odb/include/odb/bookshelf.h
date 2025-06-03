// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2019-2025, The OpenROAD Authors

#pragma once

#include "odb.h"
namespace utl {
class Logger;
}
namespace odb {

class bookshelfout_impl;
class dbBlock;

class bookshelfout
{
  bookshelfout_impl* _writer;

 public:
  bookshelfout(utl::Logger* logger);
  ~bookshelfout();

  bool writeBlock(dbBlock* block, const char* filename);
};

}  // namespace odb
