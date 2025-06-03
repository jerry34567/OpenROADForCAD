// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2019-2025, The OpenROAD Authors

#include "odb/bookshelf.h"
#include "bookshelfout_impl.h"

#include <cstdio>

#include "odb/db.h"

namespace odb {

bookshelfout::bookshelfout(utl::Logger* logger)
{
  _writer = new bookshelfout_impl(logger);
  assert(_writer);
}

bookshelfout::~bookshelfout()
{
  delete _writer;
}

bool bookshelfout::writeBlock(dbBlock* block, const char* filename)
{
  return _writer->writeBlock(block, filename);
}

}  // namespace odb
