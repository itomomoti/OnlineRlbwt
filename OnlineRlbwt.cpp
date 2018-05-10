/*!
 * Copyright (c) 2017 Tomohiro I
 *
 * This program is released under the MIT License.
 * http://opensource.org/licenses/mit-license.php
 */
/*!
 * @file OnlineRlbwt.cpp
 * @brief Online RLBWT construction.
 * @author Tomohiro I
 * @date 2018-01-27
 */
#include <stdint.h>

#include <time.h>
#include <iostream>
#include <iomanip>
#include <chrono>

#include "cmdline.h"
#include "OnlineRlbwt.hpp"
#include "DynRleForRlbwt.hpp"


using namespace itmmti;

int main(int argc, char *argv[])
{
  cmdline::parser parser;
  parser.add<std::string>("input",'i', "input file name", true);
  parser.add<std::string>("output",'o', "output file name", false);
  parser.add<bool>("check", 0, "check correctness", false, 0);
  parser.add<bool>("verbose", 'v', "verbose", false, 0);
  parser.add("help", 0, "print help");

  parser.parse_check(argc, argv);
  const std::string in = parser.get<std::string>("input");
  const std::string out = parser.get<std::string>("output");
  const bool check = parser.get<bool>("check");
  const bool verbose = parser.get<bool>("verbose");

  auto t1 = std::chrono::high_resolution_clock::now();
  std::cout << "Building RLBWT ..." << std::endl;

	std::ifstream ifs(in);

  size_t j = 0;
  const size_t step = 1000000; // print status every step characters
  size_t last_step = 0;

  using BTreeNodeT = BTreeNode<32>;
  using BtmNodeMT = BtmNodeM_StepCode<BTreeNodeT, 32>;
  using BtmMInfoT = BtmMInfo_BlockVec<BtmNodeMT, 512>;
  using BtmNodeST = BtmNodeS<BTreeNodeT, uint32_t, 8>;
  using BtmSInfoT = BtmSInfo_BlockVec<BtmNodeST, 1024>;
  using DynRleT = DynRleForRlbwt<WBitsBlockVec<1024>, Samples_Null, BtmMInfoT, BtmSInfoT>;
  OnlineRlbwt<DynRleT> rlbwt(1);

  char c; // Assume that the input character fits in char.
  unsigned char uc;

  while (ifs.peek() != std::ios::traits_type::eof()) {
    ifs.get(c);
    uc = static_cast<unsigned char>(c);
    if (verbose) {
      if(j > last_step + (step - 1)){
        last_step = j;
        std::cout << " " << j << " characters processed ..." << std::endl;
      }
    }

    rlbwt.extend(uint8_t(uc));
    ++j;
  }

  ifs.close();
  {
    auto t2 = std::chrono::high_resolution_clock::now();
    double sec = std::chrono::duration_cast<std::chrono::seconds>(t2 - t1).count();
    std::cout << "RLBWT construction done. " << sec << " sec" << std::endl;
  }

  rlbwt.printStatistics(std::cout, false);

  if (!(out.empty())) {
    t1 = std::chrono::high_resolution_clock::now();
    std::cout << "Decompressing RLBWT ..." << std::endl;
    std::ofstream ofs(out, std::ios::out);
    rlbwt.invert(ofs);
    auto t2 = std::chrono::high_resolution_clock::now();
    double sec = std::chrono::duration_cast<std::chrono::seconds>(t2 - t1).count();
    std::cout << "RLBWT decompression done. " << sec << " sec" << std::endl;
  }

  if (check) { // check correctness
    t1 = std::chrono::high_resolution_clock::now();
    std::cout << "Checking RLBWT inversion ..." << std::endl;
    std::ifstream ifssss(in);
    if (rlbwt.checkDecompress(ifssss)) {
      auto t2 = std::chrono::high_resolution_clock::now();
      double sec = std::chrono::duration_cast<std::chrono::seconds>(t2 - t1).count();
      std::cout << "RLBWT decompressed correctly. " << sec << " sec" << std::endl;
    } else {
      std::cout << "RLBWT inversion failed." << std::endl;
    }
  }
}
