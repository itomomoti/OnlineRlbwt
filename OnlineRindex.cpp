/*!
 * Copyright (c) 2018 Tomohiro I
 *
 * This program is released under the MIT License.
 * http://opensource.org/licenses/mit-license.php
 */
/*!
 * @file OnlineRindex.cpp
 * @brief Online r-index, index based on Run-length encoded Burrows–Wheeler transform (RLBWT).
 * @author Tomohiro I
 * @date 2018-05-04
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

int main(int argc, char *argv[])
{
  cmdline::parser parser;
  parser.add<std::string>("input", 'i', "input file name", true);
  parser.add<bool>("check", 0, "check correctness", false, 0);
  parser.add<bool>("verbose", 'v', "verbose", false, 0);
  parser.add("help", 0, "print help");

  parser.parse_check(argc, argv);
  const std::string in = parser.get<std::string>("input");
  const bool check = parser.get<bool>("check");
  const bool verbose = parser.get<bool>("verbose");

  auto t1 = std::chrono::high_resolution_clock::now();
  std::cout << "R-index constructing..." << std::endl;

  std::ifstream ifs(in);

  const size_t step = 1000000; // Print status every step characters.
  size_t last_step = 0;

  using BTreeNodeT = BTreeNode<32>; // BTree arity = {16, 32, 64, 128}
  using BtmNodeMT = BtmNodeM_StepCode<BTreeNodeT, 32>; // BtmNode arity in {16, 32, 64, 128}.
  using BtmMInfoT = BtmMInfo_BlockVec<BtmNodeMT, 512>; // Each block has 512 btmNodeM.
  using BtmNodeST = BtmNodeS<BTreeNodeT, uint32_t, 8>; // CharT = uint32_t. BtmNode arity = {4, 8, 16, 32, 64, 128}.
  using BtmSInfoT = BtmSInfo_BlockVec<BtmNodeST, 1024>; // Each block has 1024 btmNodeS.
  using DynRleT = DynRleForRlbwt<WBitsBlockVec<1024>, Samples_WBitsBlockVec<1024>, BtmMInfoT, BtmSInfoT>;
  using BtmNodeInSucc = BtmNodeForPSumWithVal<32>; // BtmNode arity = {16, 32, 64, 128}.
  using DynSuccT = DynSuccForRindex<BTreeNodeT, BtmNodeInSucc>;
  using RindexT = OnlineRlbwtIndex<DynRleT, DynSuccT>;
  RindexT rindex(1);
  SizeT pos = 0; // Current txt-pos (0base)
  char c; // Assume that the input character fits in char.
  unsigned char uc;

  while (ifs.peek() != std::ios::traits_type::eof()) {
    ifs.get(c);
    uc = static_cast<unsigned char>(c);
    if (verbose) {
      // if (pos >= 0) {
      //   std::cerr << "loop: " << pos
      //             << ", prev = " << rindex.getPrevSamplePos()
      //             << ", next = " << rindex.getNextSamplePos()
      //             << ", insert " << (int)c << "(" << c << ")" << " at " << rindex.getEndmarkerPos() << std::endl;
      // }
      if (pos > last_step + (step - 1)) {
        last_step = pos;
        std::cout << " " << pos << " characters processed..." << std::endl;
        // {//debug
        //   rindex.printDebugInfo(std::cout);
        // }
        // rindex.printStatistics(std::cout, false);
      }
    }

    rindex.extend(uc);
    // if (verbose) {
    //   if (pos > 0) {
    //     std::cout << "Status after inserting pos = " << pos << std::endl;
    //     rindex.printDebugInfo(std::cout);
    //     // rindex.printStatistics(std::cout);
    //   }
    // }
    ++pos;
  }

  ifs.close();

  auto t2 = std::chrono::high_resolution_clock::now();
  double sec = std::chrono::duration_cast<std::chrono::seconds>(t2 - t1).count();
  std::cout << "R-index construction done. " << sec << " sec" << std::endl;
  rindex.printStatistics(std::cout, false);

  if (check) { // check correctness
    t1 = std::chrono::high_resolution_clock::now();
    std::cout << "Checking RLBWT inversion..." << std::endl;
    std::ifstream ifssss(in);
    if (rindex.checkDecompress(ifssss)) {
      auto t2 = std::chrono::high_resolution_clock::now();
      double sec = std::chrono::duration_cast<std::chrono::seconds>(t2 - t1).count();
      std::cout << "RLBWT decompressed correctly. " << sec << " sec" << std::endl;
    } else {
      std::cout << "RLBWT inversion failed." << std::endl;
    }

    std::cout << "Checking correctness of data structures..." << std::endl;
    rindex.printDebugInfo(std::cout);
    std::cout << "Done." << std::endl;
  }

  return 0;
}
