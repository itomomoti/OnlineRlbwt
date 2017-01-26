#include "DynRLE.hpp"


#define RLBWT_TEST
#ifdef RLBWT_TEST
#include <time.h>
#include <iostream>
#include <iomanip>
#include <chrono>
#include "cmdline.h"

//
// $ g++ -O3 -DNDEBUG -Wall -std=c++14 -mavx -c ../../Basics/BitsUtil.cpp
// $ g++ -O3 -DNDEBUG -Wall -std=c++14 -mavx -o DynRLE.out DynRLE.cpp ../../Basics/BitsUtil.o
// $ ./DynRLE.out -n 1000 -r 28
//
int main(int argc, char *argv[])
{
  cmdline::parser parser;
  parser.add<std::string>("input",'i', "input file name", true);
  parser.add<std::string>("output",'o', "output file name", true);
  parser.add("help", 0, "print help");

  parser.parse_check(argc, argv);
  const std::string in = parser.get<std::string>("input");
  const std::string out = parser.get<std::string>("output");

  auto t1 = std::chrono::high_resolution_clock::now();

	std::ifstream ifs(in);
	std::ofstream os(out);

  size_t j = 0;
  const size_t step = 1000000;	//print status every step characters
  long int last_step = 0;

  std::cout << "Building RLBWT ..." << std::endl;

  DynRLBWT<DynRLE<32> > rlbwt(16);

  { // construct
    char c;
    while (ifs.get(c)) {
      if(j > last_step + (step - 1)){
        last_step = j;
        std::cout << " " << j << " characters processed ..." << std::endl;
      }

      rlbwt.extend(uint8_t(c));
      ++j;
    }
  }

  // { // 
  //   for (size_t i = 0; i < rlbwt.getLenWithTerminator(); ++i) {
  //     uint64_t pos = i;
  //     uint64_t idxM = rlbwt.searchPosM(pos);
  //     char c = rlbwt.getCharFromIdxM(idxM);
  //     os.put(c);
  //   }
  // }

  { // inverse
    uint64_t pos = 0;
    for (size_t i = 0; i < rlbwt.getLenWithTerminator() - 1; ++i) {
      uint64_t ch64 = rlbwt[pos];
      os.put(char(ch64));
      pos = rlbwt.totalRank(uint8_t(ch64), pos);
    }
  }

  auto t2 = std::chrono::high_resolution_clock::now();

	ifs.close();
	os.close();

  double sec = std::chrono::duration_cast<std::chrono::seconds>(t2 - t1).count();
  std::cout << "done. " << sec << " sec" << std::endl;
  rlbwt.printStatictics(std::cout);

	size_t bitsize = rlbwt.calcMemBytes() * 8;
	std::cout << " Size of the structures (bits): " << bitsize << std::endl;
	std::cout << " Size of the structures (Bytes): " << bitsize/8 << std::endl;
	std::cout << " Size of the structures (KB): " << (bitsize/8)/1024 << std::endl;
	std::cout << " Size of the structures (MB): " << ((bitsize/8)/1024)/1024 << std::endl;
}
#endif





// #define TEST_CORRECT_RANDOM_
#ifdef TEST_CORRECT_RANDOM_
#include <time.h>
#include <iostream>
#include <iomanip>
#include <chrono>
#include "cmdline.h"

//
// $ g++ -O3 -DNDEBUG -Wall -std=c++14 -mavx -c ../../BitsUtil.cpp
// $ g++ -O3 -DNDEBUG -Wall -std=c++14 -mavx -o DynRLE.out DynRLE.cpp ../../Basics/BitsUtil.o
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

