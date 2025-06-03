// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2019-2025, The OpenROAD Authors

#include "bookshelfout_impl.h"

#include <sys/stat.h>

#include <algorithm>
#include <cstdint>
#include <cstdio>
#include <limits>
#include <map>
#include <optional>
#include <set>
#include <string>
#include <vector>

#include "odb/db.h"
#include "odb/dbMap.h"
#include "odb/bookshelf.h"
#include "utl/Logger.h"
#include "utl/ScopedTemporaryFile.h"

namespace odb {

namespace {

template <typename T>
std::vector<T*> sortedSet(dbSet<T>& to_sort)
{
  std::vector<T*> sorted(to_sort.begin(), to_sort.end());
  std::sort(sorted.begin(), sorted.end(), [](T* a, T* b) {
    return a->getName() < b->getName();
  });
  return sorted;
}

const char* defOrient(const dbOrientType& orient)
{
  switch (orient.getValue()) {
    case dbOrientType::R0:
      return "N";

    case dbOrientType::R90:
      return "W";

    case dbOrientType::R180:
      return "S";

    case dbOrientType::R270:
      return "E";

    case dbOrientType::MY:
      return "FN";

    case dbOrientType::MYR90:
      return "FE";

    case dbOrientType::MX:
      return "FS";

    case dbOrientType::MXR90:
      return "FW";
  }

  return "N";
}

}

bool bookshelfout_impl::writeBlock(dbBlock* block, const char* filename)
{
  _dist_factor
      = (double) block->getDefUnits() / (double) block->getDbUnitsPerMicron();
  utl::FileHandler fileHandler(filename);
  _out = fileHandler.getFile();

  if (_out == nullptr) {
    _logger->warn(
        utl::ODB, 273, "Cannot open Bookshelf file ({}) for writing", filename);
    return false;
  }

  fprintf(_out, "UCLA pl 1.0\n\n");

  dbSet<dbInst> insts = block->getInsts();

  // Sort the components for consistent output
  for (dbInst* inst : insts) {
    writeInst(inst);
  }

  dbSet<dbBTerm> bterms = block->getBTerms();

  for (dbBTerm* bterm : sortedSet(bterms)) {
    writeBTerm(bterm);
  }
  return true;
}

void bookshelfout_impl::writeInst(dbInst* inst)
{
  std::string iname = inst->getName();
  fprintf(_out, "%s", iname.c_str());

  int x, y;
  inst->getLocation(x, y);
  x = defdist(x);
  y = defdist(y);

  const char* orient = defOrient(inst->getOrient());
  dbPlacementStatus status = inst->getPlacementStatus();

  switch (status.getValue()) {
    case dbPlacementStatus::SUGGESTED:
    case dbPlacementStatus::PLACED: {
      fprintf(_out, " %d %d : %s", x, y, orient);
      break;
    }

    case dbPlacementStatus::LOCKED:
    case dbPlacementStatus::FIRM: {
      fprintf(_out, " + FIXED ( %d %d ) : %s", x, y, orient);
      break;
    }

    default:
      printf("Unknown placement status: %d\n", status.getValue());
      break;
  }
  fprintf(_out, "\n");
}

void bookshelfout_impl::writeBTerm(dbBTerm* bterm)
{
  dbNet* net = bterm->getNet();
  if (net && !net->isSpecial()) {
    dbSet<dbBPin> bpins = bterm->getBPins();

    if (bpins.size() != 0) {
      int cnt = 0;

      dbSet<dbBPin>::iterator itr;

      for (itr = bpins.begin(); itr != bpins.end(); ++itr) {
        writeBPin(*itr, cnt++);
      }
      return;
    }
  }
}

void bookshelfout_impl::writeBPin(dbBPin* bpin, int cnt)
{
  dbBTerm* bterm = bpin->getBTerm();
  std::string bname = bterm->getName();
  fprintf(_out, "%s", bname.c_str());

  bool isFirst = true;
  int dw, dh, x = 0, y = 0;
  int xMin, yMin, xMax, yMax;
  for (dbBox* box : bpin->getBoxes()) {
    dw = defdist(int(box->getDX() / 2));
    dh = defdist(int(box->getDY() / 2));

    if (isFirst) {
      isFirst = false;
      x = defdist(box->xMin()) + dw;
      y = defdist(box->yMin()) + dh;
    }

    xMin = defdist(box->xMin()) - x;
    yMin = defdist(box->yMin()) - y;
    xMax = defdist(box->xMax()) - x;
    yMax = defdist(box->yMax()) - y;
  }
  // if (y + yMin == 0) {
  //   y = 0;
  // }
  // else {
  //   y = y - yMin;
  // }

  dbPlacementStatus status = bpin->getPlacementStatus();

  switch (status.getValue()) {
    case dbPlacementStatus::SUGGESTED:
    case dbPlacementStatus::PLACED: {
      fprintf(_out, " %d %d : N /FIXED_NI\n", x, y);
      break;
    }

    case dbPlacementStatus::LOCKED:
    case dbPlacementStatus::FIRM: {
      fprintf(_out, " + FIXED ( %d %d ) : N /FIXED_NI\n", x, y);
      break;
    }

    default:
      printf("Unknown placement status: %d\n", status.getValue());
      break;
  }
}

}  // namespace odb
