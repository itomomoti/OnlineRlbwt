#include "DynRLE.hpp"


// #define TEST_CORRECT_RANDOM_
#ifdef TEST_CORRECT_RANDOM_
#include <time.h>
#include <iostream>
#include <iomanip>
#include <chrono>
#include "cmdline.h"

//
// $ g++ -O3 -DNDEBUG -Wall -std=c++14 -mnative -c ../../BitsUtil.cpp
// $ g++ -O3 -DNDEBUG -Wall -std=c++14 -mnative -o DynRLE.out DynRLE.cpp ../../Basics/BitsUtil.o
// $ ./DynRLE.out -n 1000 -r 28
//
int main(int argc, char *argv[])
{
  cmdline::parser parser;
  parser.add<uint64_t>("num",'n', "num of random inseart", true);
  parser.add<int>("seed",'s', "seed of random values", true);
  // parser.add<uint8_t>("width", 'w', "bit width for packed vector", false, 8, cmdline::range(1, 64));
  // parser.add<uint64_t>("jump", 'j', "jump for random access", false, 38201); // 38201 is a prime number
  // parser.add<uint64_t>("val", 'v', "value to write (should fit in w bits)", false, 1);
  // parser.add<uint8_t>("dummy", 0, "dummy argument (do not input this)", false, 0);
  parser.add("help", 0, "print help");

  parser.parse_check(argc, argv);
  const uint64_t num = parser.get<uint64_t>("num");
  const int seed = parser.get<int>("seed");
  // const uint8_t w = parser.get<uint8_t>("width");
  // const uint64_t jump = parser.get<uint64_t>("jump");
  // uint64_t val = parser.get<uint64_t>("val");
  // const uint8_t dummy = parser.get<uint8_t>("dummy");

  // const uint64_t idxMask = bits::UINTW_MAX(bits::bitSize(num));
  // std::cout << "jump: " << jump << std::endl;

  // assert(num > 2);
  srand(seed);

  DynRLE<4> drle(1);

  auto t1 = std::chrono::high_resolution_clock::now();
  {
    for (uint64_t i = 0; i < num; ++i) {
      const uint64_t totalLen = drle.getSumOfWeight();
      const uint64_t ch = rand() % 256;
      const uint64_t insPos = rand() % (totalLen+1);
      const uint64_t exponent = rand() % 32;

      drle.insertRunWithoutReturn(ch, exponent, insPos);
      if ((i+1) % 1000000 == 0) {
        debugstream1 << i+1 << " runs inserted: ";
        auto t2 = std::chrono::high_resolution_clock::now();
        double msec = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count();
        std::cout << msec << " msec" << std::endl;
        drle.printStatictics(debugstream1);
      }
    }
  }
  {
    debugstream1 << num << " runs inserted: ";
    auto t2 = std::chrono::high_resolution_clock::now();
    double msec = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count();
    std::cout << msec << " msec" << std::endl;
    drle.printStatictics(debugstream1);
  }
  {
    debugstream1 << "rank/select check..." << std::endl;
    t1 = std::chrono::high_resolution_clock::now();
    for (const auto * rootS = drle.getFstRootS(); reinterpret_cast<uintptr_t>(rootS) != NOTFOUND; rootS = drle.getNextRootS(rootS)) {
      const uint64_t chNum = rootS->getSumOfWeight();
      const uint64_t ch = drle.getCharFromNodeS(rootS);
      uint64_t prevPos = 0;
      for (uint64_t i = 1; i <= chNum; ++i) {
        auto sel = drle.select(ch, i);
        auto rank = drle.rank(ch, sel, false);
        auto rank0 = (sel) ? drle.rank(ch, sel-1, false) : 0;
        if (rank - rank0 != 1) {
          debugstream1 << "err: ch=" << ch << ", i=" << i << ", rank=" << rank << ", select=" << sel << std::endl;
        }
        if (sel > 0 && prevPos >= sel) {
          debugstream1 << "err: ch=" << ch << ", i=" << i << ", prevPos=" << prevPos << ">= select=" << sel << std::endl;
        }
        prevPos = sel;
      }
    }
    auto t2 = std::chrono::high_resolution_clock::now();
    double msec = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count();
    std::cout << "done. " << msec << " msec" << std::endl;
  }
  {
    debugstream1 << "total rank check..." << std::endl;
    t1 = std::chrono::high_resolution_clock::now();
    uint64_t inc = 1;
    for (const auto * rootS = drle.getFstRootS(); reinterpret_cast<uintptr_t>(rootS) != NOTFOUND; rootS = drle.getNextRootS(rootS)) {
      const uint64_t chNum = rootS->getSumOfWeight();
      const uint64_t ch = drle.getCharFromNodeS(rootS);
      for (uint64_t i = 1; i <= chNum; ++i) {
        const uint64_t sel = drle.select(ch, i);
        auto totalRank = drle.rank(ch, sel, true); // total rank
        if (totalRank != inc) {
          debugstream1 << "err: select(" << ch << ", " << i << ")=" << sel << ", totalRank(sel)=" << totalRank << "should be " << inc << std::endl;
        }
        ++inc;
      }
    }
    auto t2 = std::chrono::high_resolution_clock::now();
    double msec = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count();
    std::cout << "done. " << msec << " msec" << std::endl;
  }
}
#endif



// #define DEBUG_SIMPLE
#ifdef DEBUG_SIMPLE
    for (uint64_t i = 0; i < 10; ++i) {
      drle.pushbackRun(i, 2);
      //    drle.printRuns();
    }
    drle.printDebugInfo(debugstream1);
    drle.printStatictics(debugstream1);
    debugstream1 << drle.calcMemBytes() << std::endl;
    // for (int64_t i = 10; i >= 0; --i) {
    //   uint64_t pos = i * 2 + 1;
    //   drle.insertRun(50, 20, pos);
    //   //    drle.printRuns();
    // }
    for (int64_t i = 0; i < 10; ++i) {
      uint64_t pos = i * 12 + 1;
      drle.insertRun(5, 10, pos);
      //    drle.printRuns();
    }
    debugstream1 << std::endl;
    drle.printDebugInfo(debugstream1);
    drle.printStatictics(debugstream1);

#endif

