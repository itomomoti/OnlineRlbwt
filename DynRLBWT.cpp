#define RLBWT_TEST
#ifdef RLBWT_TEST
#include <time.h>
#include <iostream>
#include <iomanip>
#include <chrono>

#include "cmdline.h"
#include "DynRLE.hpp"
#include "DynRLBWT.hpp"


#define MYDEBUG
#ifdef MYDEBUG
#define debugstream1 std::cout
#define debugstream2 std::cerr
#endif

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
    for (size_t i = 0; i < rlbwt.getLenWithEm() - 1; ++i) {
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
