/*!
 * Copyright (c) 2017 Tomohiro I
 *
 * This program is released under the MIT License.
 * http://opensource.org/licenses/mit-license.php
 */
/*!
 * @file OnlineLz77ViaRlbwt.cpp
 * @brief Online LZ77 computation via online RLBWT.
 * @author Tomohiro I
 * @date 2017-10-12
 */
#include <stdint.h>

#include <time.h>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <chrono>

#include "cmdline.h"
#include "OnlineLz77ViaRlbwt.hpp"
#include "DynRleForRlbwt.hpp"


using namespace itmmti;
using SizeT = uint32_t; // Text length should fit in SizeT.

//
// $ ./LZ77ViaOnlineRlbwt -i inputfilename -o outputfilename
//
int main(int argc, char *argv[])
{
  cmdline::parser parser;
  parser.add<std::string>("input", 'i', "input file name", true);
  parser.add<std::string>("output", 'o', "output file name", true);
  parser.add<bool>("verbose", 'v', "verbose", false, 0);
  parser.add("help", 0, "print help");

  parser.parse_check(argc, argv);
  const std::string in = parser.get<std::string>("input");
  const std::string out = parser.get<std::string>("output");
  const bool verbose = parser.get<bool>("verbose");

  auto t1 = std::chrono::high_resolution_clock::now();

  std::ifstream ifs(in);
  std::ofstream ofs(out);

  const size_t step = 1000000; // Print status every step characters.
  size_t last_step = 0;
  std::cout << "LZ77 Parsing ..." << std::endl;

  OnlineLz77ViaRlbwt<DynRleForRlbwt<32, 32> > rlbwt(128);
  SizeT pos = 0; // Current txt-pos (0base)
  SizeT l = 0; // Length of current LZ phrase prefix
  SizeT z = 0; // LZ phrase counter
  bwttracker tracker = {0, 1, 0}; // BWT tracker for current LZ phrase prefix.
  char c; // Assume that the input character fits in char.
  unsigned char uc;

  // {//test
  //   if (reinterpret_cast<void *>(0) == nullptr) {
  //     std::cout << "a" << std::endl;
  //   } else {
  //     std::cout << "b" << std::endl;
  //   }
  //   return 0;
  // }

  while (ifs.peek() != std::ios::traits_type::eof()) {
    ifs.get(c);
    uc = static_cast<unsigned char>(c);
    if (verbose) {
      // std::cout << pos << ": z = " << z << ", l = " << l << ", ref = " << std::get<2>(tracker) - l << ", succ = " << rlbwt.getSuccSamplePos()
      //           << ", tracker = (" << std::get<0>(tracker) << ", " << std::get<1>(tracker) << ", " << std::get<2>(tracker) << "), "
      //           << "insert " << (int)c << "(" << c << ")" << " at " << rlbwt.getEndmarkerPos() << std::endl;
      if (pos > last_step + (step - 1)) {
        last_step = pos;
        std::cout << " " << pos << " characters processed ..." << std::endl;
        // {//debug check if text can be decompressed correctly from rlbwt
        //   std::ifstream ifssss(in);
        //   rlbwt.checkDecompress(ifssss);
        // }
      }
    }

    if (rlbwt.lfMap(tracker, uc)) {
      ++l;
      ++(std::get<1>(tracker)); // New suffix falls inside range
    } else { // Extension failed: End of an LZ factor.
      // Output an LZ factor.
      ++z;
      SizeT beg = static_cast<SizeT>(std::get<2>(tracker) - l);
      ofs.write(reinterpret_cast<char *>(&beg), sizeof(SizeT));
      ofs.write(reinterpret_cast<char *>(&l), sizeof(SizeT));
      ofs.write(&c, 1);
      // if (verbose) {
      //   std::cout << "LZ[" << z << "] = (" << beg << ", " << l << ", " << c << ")" << std::endl;
      // }
      // Reset variables.
      l = 0;
      tracker = {0, rlbwt.getLenWithEndmarker() + 1, 0}; // +1 because we have not inserted "ch".
    }

    rlbwt.extend(uc);
    ++pos;
    if (std::get<0>(tracker) == rlbwt.getEndmarkerPos()) {
      std::get<2>(tracker) = rlbwt.getSuccSamplePos();
    }

    if (verbose) {
      // if (pos >= 677) {
      //   std::cout << "Status after inserting pos = " << pos - 1 << std::endl;
      //   rlbwt.printDebugInfo(std::cout);
      //   {//debug check if text can be decompressed correctly from rlbwt
      //     std::ifstream ifssss(in);
      //     rlbwt.checkDecompress(ifssss);
      //   }
      // }
      // if (pos > 680) {
      //   exit(1);
      // }
    }
  }
  if (l) {
    ++z;
    SizeT beg = static_cast<SizeT>(std::get<2>(tracker) - l);
    SizeT len = l - 1;
    ofs.write(reinterpret_cast<char *>(&beg), sizeof(SizeT));
    ofs.write(reinterpret_cast<char *>(&len), sizeof(SizeT));
    ofs.write(&c, 1);
    // if (verbose) {
    //   std::cout << "LZ[" << z << "] = (" << beg << ", " << l-1 << ", " << c << ")" << std::endl;
    // }
  }

  ifs.close();
  ofs.close();

  auto t2 = std::chrono::high_resolution_clock::now();
  double sec = std::chrono::duration_cast<std::chrono::seconds>(t2 - t1).count();
  std::cout << "LZ compression done. " << sec << " sec" << std::endl;
  std::cout << "Number of factors z = " << z << std::endl;
  rlbwt.printStatictics(std::cout);

  return 0;
}
