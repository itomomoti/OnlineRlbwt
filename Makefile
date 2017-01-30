CXX = g++
CXXFLAGS = -std=c++14 -march=native -O3 -DNDEBUG -W -Wall -Wno-deprecated
CXXFLAGS += -MD
LINKFLAGS = -lm

BASICS_DIR = ../../Basics/

SRCS =\
	DynRLE.cpp\
	$(BASICS_DIR)BitsUtil.cpp\

SRCS1 = DynRLBWT.cpp\

OBJS = $(SRCS:%.cpp=%.o)

OBJS1 = $(SRCS1:%.cpp=%.o)

all: DynRLBWT

DynRLBWT:  $(OBJS) $(OBJS1)
	$(CXX) $(CXXFLAGS) $(OTHERFLAGS) $(OBJS) $(OBJS1) $(LINKFLAGS) -o DynRLBWT.out

debug:
	make all CXXFLAGS="-std=c++14 -march=native -g3 -W -Wall -MD"

clean:
	rm -f *.o *~
	rm -f *.d

-include *.d
