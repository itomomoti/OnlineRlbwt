/*!
 * Copyright (c) 2017 Tomohiro I
 *
 * This program is released under the MIT License.
 * http://opensource.org/licenses/mit-license.php
 */
/*!
 * @file DecompressLz77.cpp
 * @brief Decompress LZ77.
 * @author Tomohiro I
 * @date 2017-10-12
 */
#include <stdint.h>

#include <time.h>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <chrono>
#include <vector>

#include "cmdline.h"


using SizeT = uint32_t; // Text length should fit in SizeT.
using txt = std::vector<char>;

int main(int argc, char *argv[])
{
  cmdline::parser parser;
  parser.add<std::string>("input",'i', "input file name", true);
  parser.add<std::string>("output",'o', "output file name", true);
  parser.add<bool>("verbose",'v', "verbose", false, 0);
  parser.add("help", 0, "print help");

  parser.parse_check(argc, argv);
  const std::string in = parser.get<std::string>("input");
  const std::string out = parser.get<std::string>("output");
  const bool verbose = parser.get<bool>("verbose");

  auto t1 = std::chrono::high_resolution_clock::now();

  std::cout << "LZ77 Decompressing ..." << std::endl;
  txt T;

  SizeT n = 0; // Length of text.
  SizeT z = 0; // LZ phrase counter
  {
    std::ifstream ifs(in, std::ios::in | std::ios::binary);
    SizeT beg;
    SizeT len;
    char c;
    while (ifs.peek() != std::ios::traits_type::eof()) {
      ifs.read((char *) &beg, sizeof(SizeT));
      ifs.read((char *) &len, sizeof(SizeT));
      ifs.read((char *) &c, 1);

      n += len + 1;
      ++z;

      // if (verbose) {
      //   std::cout << "LZ[" << z << "] = (" << beg << ", " << len << ", " << c << ")" << std::endl;
      // }
      for (SizeT i = 0; i < len; ++i) {
        T.push_back(T[beg + i]);
      }
      T.push_back(c);
    }
    ifs.close();
  }
  {
    std::ofstream ofs(out);
    for (SizeT i = 0; i < T.size(); ++i) {
      ofs.write(&(T[i]), 1);
    }
    ofs.close();
  }

  auto t2 = std::chrono::high_resolution_clock::now();
  double sec = std::chrono::duration_cast<std::chrono::seconds>(t2 - t1).count();
  std::cout << "LZ77 Decompression done. " << sec << " sec" << std::endl;
  std::cout << "Text length n = " << n << std::endl;
  std::cout << "Number of factors z = " << z << std::endl;
}
