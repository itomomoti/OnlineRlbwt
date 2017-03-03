#define RLBWT_TEST
#ifdef RLBWT_TEST
#include <time.h>
#include <iostream>
#include <iomanip>
#include <chrono>

#include "cmdline.h"
#include "DynRLE.hpp"
#include "OnlineRLBWT.hpp"


//
// $ ./DynRLBWT.out -i inputfilename -o outputfilename
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
	std::ofstream ofs(out);

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
    auto t2 = std::chrono::high_resolution_clock::now();
    double sec = std::chrono::duration_cast<std::chrono::seconds>(t2 - t1).count();
    std::cout << "construct done. " << sec << " sec" << std::endl;
    rlbwt.printStatictics(std::cout);
  }


  // { // 
  //   for (size_t i = 0; i < rlbwt.getLenWithTerminator(); ++i) {
  //     uint64_t pos = i;
  //     uint64_t idxM = rlbwt.searchPosM(pos);
  //     char c = rlbwt.getCharFromIdxM(idxM);
  //     ofs.put(c);
  //   }
  // }

  { // inverse
    rlbwt.invert(ofs);
  }

  auto t2 = std::chrono::high_resolution_clock::now();

	ifs.close();
	ofs.close();

  double sec = std::chrono::duration_cast<std::chrono::seconds>(t2 - t1).count();
  std::cout << "construct and string inversion done. " << sec << " sec" << std::endl;
  rlbwt.printStatictics(std::cout);

	size_t bitsize = rlbwt.calcMemBytes() * 8;
	std::cout << " Size of the structures (bits): " << bitsize << std::endl;
	std::cout << " Size of the structures (Bytes): " << bitsize/8 << std::endl;
	std::cout << " Size of the structures (KB): " << (bitsize/8)/1024 << std::endl;
	std::cout << " Size of the structures (MB): " << ((bitsize/8)/1024)/1024 << std::endl;
}
#endif
