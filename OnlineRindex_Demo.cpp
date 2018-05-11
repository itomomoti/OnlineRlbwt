/*!
 * Copyright (c) 2018 Tomohiro I
 *
 * This program is released under the MIT License.
 * http://opensource.org/licenses/mit-license.php
 */
/*!
 * @file OnlineRindex_Demo.cpp
 * @brief Demonstration for online r-index, index based on Run-length encoded Burrows–Wheeler transform (RLBWT).
 * @author Tomohiro I
 * @date 2018-05-10
 */
#include <stdint.h>

#include <time.h>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <chrono>

#include "cmdline.h"
#include "OnlineRindex.hpp"
#include "DynRleForRlbwt.hpp"
#include "DynSuccForRindex.hpp"


using namespace itmmti;
using SizeT = uint64_t; // Text length should fit in SizeT.

template<typename RindexT>
void searchOnRindex(const RindexT & rindex, std::string message)
{
  //// Pattern search Demo
  const uint64_t lenWithoutEm = rindex.getLenWithEndmarker() - 1;
  while (true) {
    constexpr static uint64_t bufsize = 512;
    char buffer[bufsize] = "";
    std::cout << message << std::endl;
    if (std::fgets(buffer, bufsize, stdin) == NULL || buffer[0] == '\n') {
      break;
    }
    uint64_t len = std::strlen(buffer);
    if (len > 0) {
      if (buffer[len - 1] == '\n') {
        buffer[--len] = '\0';
      } else {
        while (getchar() != '\n') {};
      }
    }
    std::cout << "\"" << buffer << "\"" << std::endl;

    //// Counting...
    auto t1 = std::chrono::high_resolution_clock::now();
    auto tracker = rindex.getInitialPatTracker();
    bool match = true;
    uint64_t plen = 0;
    for (unsigned char c = buffer[plen]; c != '\0'; c = buffer[++plen]) {
      match = rindex.lfMap(tracker, c);
      if (!match) break;
    }
    auto t2 = std::chrono::high_resolution_clock::now();
    double microsec = std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count();
    std::cout << "Counting: done in " << microsec << " micro sec. Len = " << len << ". Prefix match length = " << plen << ". ";

    if (match) {
      const auto numOcc = rindex.getNumOcc(tracker);
      std::cout << "NumOcc = " << numOcc << ". BwtInterval = [" << std::get<0>(tracker) << ".."
                << std::get<1>(tracker) << ")" << std::endl;
      //// Note that the returned occ is end position (exclusive) of pattern.
      //// To get beginning position, subtract pattern length.
      unsigned char c;
      do {
        std::cout << "Print them all (y/n)? Or locate without printing (l): ";
        c = getchar();
        while (getchar() != '\n') {};
      } while (!(c == 'y' || c == 'n' || c == 'l'));
      //// Locating...
      if (c == 'l') { // Locating without printing.
        t1 = std::chrono::high_resolution_clock::now();
        auto endPos = rindex.calcFstOcc(tracker);
        for (auto num = numOcc; num; --num) {
          endPos = rindex.calcNextPos(endPos);
        }
        t2 = std::chrono::high_resolution_clock::now();
        microsec = std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count();
        std::cout << "Locating: done in " << microsec << " micro sec. "
                  << microsec / numOcc << " micro sec. each." << std::endl;
        if (!rindex.isReady()) { // dummy code to prevent optimization delete the locating codes
          std::cout << endPos << std::endl;
        }
      } else if (c == 'y') { // Locating with printing.
        t1 = std::chrono::high_resolution_clock::now();
        auto endPos = rindex.calcFstOcc(tracker);
        for (auto num = numOcc; num; --num) {
          std::cout << endPos - len << ", ";
          endPos = rindex.calcNextPos(endPos);
        }
        std::cout << std::endl;
        t2 = std::chrono::high_resolution_clock::now();
        microsec = std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count();
        std::cout << "Locating: done in " << microsec << " micro sec. " << microsec / numOcc << " micro sec. each." << std::endl;
      }
    } else {
      std::cout << "Pttern does not match." << std::endl;
    }
  }
}


int main(int argc, char *argv[])
{
  cmdline::parser parser;
  parser.add<std::string>("input", 'i', "input file name", true);
  parser.add<size_t>("step", 's', "number of characters to index in a single step", false, 1000000);
  parser.add<bool>("verbose", 'v', "verbose", false, 0);
  parser.add("help", 0, "print help");

  parser.parse_check(argc, argv);
  const std::string in = parser.get<std::string>("input");
  const bool verbose = parser.get<bool>("verbose");
  const size_t step = parser.get<size_t>("step");

  std::cout << "R-index constructing..." << std::endl;

  std::ifstream ifs(in);

  size_t last_step = 0;

  using BTreeNodeT = BTreeNode<16>; // BTree arity = {16, 32, 64, 128}
  using BtmNodeMT = BtmNodeM_StepCode<BTreeNodeT, 32>; // BtmNode arity in {16, 32, 64, 128}.
  using BtmMInfoT = BtmMInfo_BlockVec<BtmNodeMT, 512>; // Each block has 512 btmNodeM.
  using BtmNodeST = BtmNodeS<BTreeNodeT, uint32_t, 8>; // CharT = uint32_t. BtmNode arity = {4, 8, 16, 32, 64, 128}.
  using BtmSInfoT = BtmSInfo_BlockVec<BtmNodeST, 1024>; // Each block has 1024 btmNodeS.
  using DynRleT = DynRleForRlbwt<WBitsBlockVec<1024>, Samples_WBitsBlockVec<1024>, BtmMInfoT, BtmSInfoT>;
  using BtmNodeInSucc = BtmNodeForPSumWithVal<16>; // BtmNode arity = {16, 32, 64, 128}.
  using DynSuccT = DynSuccForRindex<BTreeNodeT, BtmNodeInSucc>;
  using RindexT = OnlineRlbwtIndex<DynRleT, DynSuccT>;
  RindexT rindex(1);
  SizeT pos = 0; // Current txt-pos (0base)
  char c; // Assume that the input character fits in char.
  unsigned char uc;

  while (ifs.peek() != std::ios::traits_type::eof()) {
    ifs.get(c);
    uc = static_cast<unsigned char>(c);

    if (pos > last_step + (step - 1)) {
      if (verbose) {
        rindex.printStatistics(std::cout, false);
      }
      last_step = pos;
      const size_t totalBytes = rindex.calcMemBytes(true);
      std::cout << " " << pos << " characters indexed in "
                << totalBytes << " bytes = "
                << (double)(totalBytes) / 1024 << " KiB = "
                << ((double)(totalBytes) / 1024) / 1024 << " MiB." << std::endl;
      searchOnRindex(rindex, "Type a pattern to search. Or enter empty string to continue indexing.");
      std::cout << "Quitted searching phase and continue indexing next " << step << " characters..." << std::endl;
    }

    rindex.extend(uc);
    ++pos;
  }

  ifs.close();

  const size_t totalBytes = rindex.calcMemBytes(true);
  std::cout << " " << pos << " characters indexed in "
            << totalBytes << " bytes = "
            << (double)(totalBytes) / 1024 << " KiB = "
            << ((double)(totalBytes) / 1024) / 1024 << " MiB." << std::endl;
  searchOnRindex(rindex, "Type a pattern to search. Or enter empty string to quit.");
  std::cout << "Quitted." << std::endl;
  if (verbose) {
    rindex.printStatistics(std::cout, false);
  }

  return 0;
}
