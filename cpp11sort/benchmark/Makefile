CXX = g++
#CXX = icpc 

CXXFLAGS += -std=c++14
CXXFLAGS += -O3 -DNDEBUG
CXXFLAGS += -I../include

LDFLAGS += -lpthread

# GNU libstc++ parallel mode
CXXFLAGS += -DWITH_STDCPP_PM
CXXFLAGS += -fopenmp
LDFLAGS += -fopenmp

# Intel TBB 
CXXFLAGS += -DWITH_TBB
CXXFLAGS += -I$(TBBROOT)/include
LDFLAGS += -ltbb
LDFLAGS += -L$(TBBROOT)/lib/intel64/gcc4.7

# AQsort
#CXXFLAGS += -DWITH_AQSORT
#CXXFLAGS += -I$(AQSORT_ROOT)/include
#CXXFLAGS += -fopenmp
#LDFLAGS += -fopenmp

# Intel Parallel STL
#CXXFLAGS += -DWITH_INTEL_PSTL

.PHONY: all
all: main

main: main.o
	$(CXX) -o $@ $< $(LDFLAGS)

main.o: main.cpp
	$(CXX) $(CXXFLAGS) -c -o $@ $<

.PHONY: clean
clean:
	rm -rf main.o main
